
## To create a fortran dynamic library from a fortran subroutine for use in regent

1. Copy the subroutine to a new file, "[subroutine-name].f90". 

2. Add `bind(c, name="[subroutine-name]_")` to the subroutine header

        subroutine [subroutine-name]([arguments]) bind(c, name="[subroutine-name]_")

3. At the start of the file, add the following code: 

        MODULE module
        [add any parameters needed by the subroutine here]
        CONTAINS

4. At the end of the file, add the following code: 

        END MODULE module

## To compile the fortran function to a dynamic library file 

1. Make sure you do compilation on a compute node 

        ssh g0003

2. Compile the object files, then compile the object files together to a dynamic library file. 

        gfortran -c [subroutine-name].f90
        gfortran -shared -o lib[subroutine-name].so [subroutine-name].o

## To use the fortran function within regent 

1. Add this function: 

       function raw_ptr_factory(ty)
            local struct raw_ptr
            {
                ptr : &ty,
                offset : int,
            }
            return raw_ptr
        end

2. Add these constructs for the types you need: 

        local raw_ptr_[type] = raw_ptr_factory([type])
        terra get_raw_ptr_[type](pr : c.legion_physical_region_t,
                                 fld : c.legion_field_id_t)
            var fa = c.legion_physical_region_get_field_accessor_array_1d(pr, fld)
            var rect : c.legion_rect_1d_t
            var subrect : c.legion_rect_1d_t
            var offsets : c.legion_byte_offset_t[2]
            var ptr = c.legion_accessor_array_1d_raw_rect_ptr(fa, rect, &subrect, offsets)
            c.legion_accessor_array_1d_destroy(fa)
            return raw_ptr_[type] { ptr = [&[type]](ptr), offset = offsets[1].offset / sizeof([type]) }
        end

2. Link the library in regent with

        regentlib.linklibrary(".../mpas-regent/fortran/lib[subroutine-name].so")
        local fortranmodule = terralib.includecstring [[
            extern void [subroutine-name]_([arguments]);
        ]]

    Note: All Fortran variables are passed by reference, so where a Fortran argument is an integer, the argument in regent is an int *.

3. Create terra and regent functions to call the fortran function. 

        terra [subroutine-name]_terra(pr_[field] : c.legion_physical_region_t,
                                      fld_[field] : c.legion_field_id_t)
            var raw[field] = get_raw_ptr_[type](0, 0, 0, pr_[field], fld_[field])
            fortranmodule.[subroutine-name]_(raw[field].ptr)
        end

        task [subroutine-name]([region] : region(ispace(int1d), example_fs))
        where
            reads writes ([region].[field])
        do
            [subroutine-name]_terra(__physical([routine].[field])[0], 
                                    __fields([routine].[field])[0])
        end

    The pointer that is passed to fortran behaves and can be manipulated similar to a pointer in c. Pointer arithmatic can be applied.

### See full example in examples.f90 and examples.rg.
Another example of Regent calling Fortran code: https://gitlab.com/StanfordLegion/legion/-/blob/master/language/examples/cholesky.rg
