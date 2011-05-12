!-------------------------------------------------------!
!----------------- Stingray motion ---------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  type(set) :: body
  open(7,file='inp.geom')
!
  write(7,'(a4)') 'body'
  body = .set.'sim_body_part1.IGS'.or.&
       (.set.'body2_part2.IGS'.and..set.plane(1,(/0, 1,0/),(/0., .633,0./),0,0)).or.&
       (.set.'body2_part3.IGS'.and..set.plane(1,(/0,-1,0/),(/0.,-.633,0./),0,0))
  body = body.map.init_raji(amp=0.2,f=1.,k=2.,peak=90.,tail=0.)
  body = body.or.(.set.'body2_tail.IGS'.and..set.plane(1,(/1,0,0/),(/3,0,0/),0,0))
  call set_write(7,body)

end program design_geom
