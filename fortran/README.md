
To create a fortran dynamic library from a fortran subroutine for use in regent
-------------------------------------------------------------------------------

Copy the subroutine to a new file, "[subroutine-name].f90". 
Add `bind(c, name="[subroutine-name]_")` to the subroutine header

    subroutine [subroutine-name]([arguments]) bind(c, name="[subroutine-name]_")

At the start of the file, add the following code: 

    MODULE module
    [add any parameters needed by the subroutine here]
    CONTAINS

At the end of the file, add the following code: 

    END MODULE module


To compile the fortran function to a dynamic library file 
---------------------------------------------------------

Make sure you do compilation on a compute node 

    ssh g0003

Compile the object files, then compile the object files together to a dynamic library file. 

    gfortran -c [subroutine-name].f90
    gfortran -shared -o lib[subroutine-name].so [subroutine-name].o

To use the fortran function within regent 
-----------------------------------------

Link the library in regent with

    local fortranmodule = terralib.includecstring [[
        extern void [subroutine-name]_([arguments]);
    ]]
    regentlib.linklibrary(".../mpas-regent/fortran/lib[subroutine-name].so")

Add these constructs to your regent code for the types you need: 

    function raw_ptr_factory(ty)
        local struct raw_ptr
        {
            ptr : &ty,
            offset : int,
        }
        return raw_ptr
    end

    local raw_ptr_int = raw_ptr_factory(int)

    terra get_raw_ptr_int(y : int, x : int, bn : int,
                          pr : c.legion_physical_region_t,
                          fld : c.legion_field_id_t)
        var fa = c.legion_physical_region_get_field_accessor_array_1d(pr, fld)
        var rect : c.legion_rect_1d_t
        var subrect : c.legion_rect_1d_t
        var offsets : c.legion_byte_offset_t[2]
        rect.lo.x[0] = y * bn
        rect.lo.x[1] = x * bn
        rect.hi.x[0] = (y + 1) * bn - 1
        rect.hi.x[1] = (x + 1) * bn - 1
        var ptr = c.legion_accessor_array_1d_raw_rect_ptr(fa, rect, &subrect, offsets)
        c.legion_accessor_array_1d_destroy(fa)
        return raw_ptr_int { ptr = [&int](ptr), offset = offsets[1].offset / sizeof(int) }
    end

Example: accessing a single value in a field in a region: 

    fspace example_fs {
        A : int,
    }

    terra subroutineA_terra(pr_A : c.legion_physical_region_t,
                            fld_A : c.legion_field_id_t)
        var rawA = get_raw_ptr_int(0, 0, 0, pr_A, fld_A)
        fortranmodule.subroutineA(rawA)
    end

    task subroutineA(example_r : region(ispace(int1d), example_fs))
    where
        reads writes (example_r.A)
    do
        subroutineA_terra(__physical(phys_tbls.A)[0], 
                          __fields(phys_tbls.A)[0])
    end

Example: accessing an array value in a field in a region 


Example of Regent calling Fortran code: https://gitlab.com/StanfordLegion/legion/-/blob/master/language/examples/cholesky.rg
