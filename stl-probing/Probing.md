# Lotus Probing Fields

#### (See Lotus/solver/oop/body.f90)

The updated body type allows to probe surface values of scalr fields on the body if this one has been defined based on an `.stl` file. To build it, we import the body module as usual

```fortran
use bodyMod, only: body
```
and then create a geometry based on this body type

```fortran
type(body) :: geom
```

## The standard way (no probing)

The standard way to define a body based on a `.stl` file is to write a function that reads-in the deisred file, and build the geometry based on this. This method will work if you do not require the `.stl` nodes (actually facet centroids) to be stored to be used later when probing fields.
```fortran
type(set) function sphere(R)

    use geom_shape

    real,intent(in) :: R
    type(model_info) :: mod_info

    mod_info%file = '../sphere.stl'     ! file to use
    mod_info%s = R*(/1.,1.,1./)         ! scale it
    mod_info%x = (/0.,0.,0./)           ! define origin
    mod_info%r = (/0.,0.,0./)           ! rotation

    sphere = model_init(mod_info)

end function sphere
```
This function return a `type(set)` that is then assigned to the geometry onto which we can add a mapping
```fortran
type(body) :: geom

geom = sphere(0.5*D).map.init_rigid(2,y)
```

## The other way
The way you have to initialize the geometry if you want to be able to probe nodal values of the flow on the body is as follows.
```fortran
subroutine sphere(geom)

    use geom_shape

    class(body),intent(inout) :: geom
    class(model),pointer :: sphere_model
    type(stlth),pointer :: ptr
    type(set) :: sets
    type(model_info) :: mod_info

    mmod_info%file = '../sphere.stl'     ! file to use
    mod_info%s = R*(/1.,1.,1./)          ! scale it
    mod_info%x = (/0.,0.,0./)            ! define origin
    mod_info%r = (/0.,0.,0./)            ! rotation

    call model_init_ptr(mod_info, sphere_model, ptr)
    sets = sphere_model
    geom = sets.map.init_rigid(2,y)      ! apply mapping
    geom%srf => ptr                      ! store nodes

end subroutine sphere
```
We now have to use a __subroutine__ (`model_init_ptr`) to initialize the geometry because in addition to returning a `model` we also need it to return a pointer to the stl facets (the `ptr` in this case). We can then assign this pointer to the body and apply mapping on the body. The geom is initialzed in the main part of the code as
```fortran
type(body) :: geom

call sphere(geom)
```

The rest of the fluid `lotus.f90` file is the same, except that now you can output nodal values of, say the pressure, every N time step
```fortran
do while(flow%time<finish)

    call geom%update(flow%time+flow%dt)   ! update position
    call flow%update(geom)

    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt, 2.*geom%pforce(flow%pressure)/(pi*D**2/4.)
    
    if(mod(flow%time,wtime)<flow%dt) then
      if(root) print '(f6.1,",",f6.3)',flow%time/D,flow%dt

      call geom%writePoints(flow%pressure,flow%time)        ! write pressure point to "surf.vtp"
      call flow%write(geom)
      call display(flow%velocity%vorticity_Z(),'out_vort',lim=20./D)
    end if

end do
```
The main loop using the standard body type would be exactly the same, without the `call geom%writePoints(...)` call. The `geom%writePoints(...)` subroutine can also write to `.csv` files instead of the standard `vtp`. To do so just add the falg `csv=.true.` in the subroutine call (this only writes the specified scalar field, and not the normals at in the `vtp`). 
