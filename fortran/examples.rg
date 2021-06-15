import "regent"
local c = regentlib.c
local format = require("std/format")

-- Example of two fortran functions imported into regent. 
-- The first function sets a single integer value.
-- The second function sets a 8-byte floating value array of size 3. 

fspace example_fs {
    a : int,
    b : double[3],
}

function raw_ptr_factory(ty)
    local struct raw_ptr
    {
        ptr : &ty,
        offset : int,
    }
    return raw_ptr
end

local raw_ptr_int = raw_ptr_factory(int)
terra get_raw_ptr_int(pr : c.legion_physical_region_t,
                      fld : c.legion_field_id_t)
    var fa = c.legion_physical_region_get_field_accessor_array_1d(pr, fld)
    var rect : c.legion_rect_1d_t
    var subrect : c.legion_rect_1d_t
    var offsets : c.legion_byte_offset_t[2]
    var ptr = c.legion_accessor_array_1d_raw_rect_ptr(fa, rect, &subrect, offsets)
    c.legion_accessor_array_1d_destroy(fa)
    return raw_ptr_int { ptr = [&int](ptr), offset = offsets[1].offset / sizeof(int) }
end

local raw_ptr_double = raw_ptr_factory(double)
terra get_raw_ptr_double(pr : c.legion_physical_region_t,
                         fld : c.legion_field_id_t)
    var fa = c.legion_physical_region_get_field_accessor_array_1d(pr, fld)
    var rect : c.legion_rect_1d_t
    var subrect : c.legion_rect_1d_t
    var offsets : c.legion_byte_offset_t[2]
    var ptr = c.legion_accessor_array_1d_raw_rect_ptr(fa, rect, &subrect, offsets)
    c.legion_accessor_array_1d_destroy(fa)
    return raw_ptr_double { ptr = [&double](ptr), offset = offsets[1].offset / sizeof(double) }
end

regentlib.linklibrary("/home/zengcs/mpas/mpas-regent/fortran/libexamples.so")

local fortranmodule = terralib.includecstring [[
    extern void subroutineA_(int *a);
    extern void subroutineB_(double *b);
]]

terra subroutineA_terra(pr_A : c.legion_physical_region_t,
                        fld_A : c.legion_field_id_t,
                        i : int)
    var rawA = get_raw_ptr_int(pr_A, fld_A)
    fortranmodule.subroutineA_(rawA.ptr + i)
end

task subroutineA(example_r : region(ispace(int1d), example_fs), 
                 i : int)
where
    reads writes (example_r.a)
do
    subroutineA_terra(__physical(example_r.a)[0], 
                      __fields(example_r.a)[0],
                      i)
end

terra subroutineB_terra(pr_B : c.legion_physical_region_t,
                        fld_B : c.legion_field_id_t,
                        i : int)
    var rawB = get_raw_ptr_double(pr_B, fld_B)
    fortranmodule.subroutineB_(rawB.ptr + i*3)
end

task subroutineB(example_r : region(ispace(int1d), example_fs), 
                 i : int)
where
    reads writes (example_r.b)
do
    subroutineB_terra(__physical(example_r.b)[0], 
                      __fields(example_r.b)[0],
                      i)
end

task main()
    var example_ispace = ispace(int1d, 5)
    var example_r = region(example_ispace, example_fs)

    -- Run subroutineA for all 5 instances of a in the region
    format.println("a: {} {} {} {} {}", example_r[0].a,
                    example_r[1].a, example_r[2].a,
                    example_r[3].a, example_r[4].a)
    for i = 0, 5 do
        subroutineA(example_r, i)
        format.println("a: {} {} {} {} {}", example_r[0].a,
                        example_r[1].a, example_r[2].a,
                        example_r[3].a, example_r[4].a)
    end

    -- Run subroutineB for both arrays of b in the region
    format.println("b: [{}, {}, {}], [{}, {}, {}]",
                    example_r[0].b[0], example_r[0].b[1],
                    example_r[0].b[2], example_r[1].b[0],
                    example_r[1].b[1], example_r[1].b[2])
    for i = 0, 2 do
        subroutineB(example_r, i)
        format.println("b: [{}, {}, {}], [{}, {}, {}]",
                        example_r[0].b[0], example_r[0].b[1],
                        example_r[0].b[2], example_r[1].b[0],
                        example_r[1].b[1], example_r[1].b[2])
    end
end
regentlib.start(main)
