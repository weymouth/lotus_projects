#!MC 1120
# Created by Tecplot 360 build 11.3.29.563
$!EXPORTSETUP EXPORTFORMAT = AVI
$!EXPORTSETUP IMAGEWIDTH = 1920
$!EXPORTSETUP ANIMATIONSPEED = 12
$!EXPORTSETUP EXPORTFNAME = 'HD.avi'

$!ADDONCOMMAND ADDONID='Extend Time MCR' 
  COMMAND='QUERY.NUMTIMESTEPS NUMTIMESTEPS'

$!LOOP |NUMTIMESTEPS|
   $!Varset |tstep| = |LOOP|
   $!LOOP |NUMFRAMES|
      $!FRAMECONTROL PUSHTOP
      $!EXTENDEDCOMMAND 
        COMMANDPROCESSORID='extend time mcr' 
        COMMAND='SET.CURTIMESTEP |tstep|'
   $!ENDLOOP
#   $!ROTATE3DVIEW THETA
#     ANGLE = 1
#     ROTATEORIGINLOCATION = DEFINEDORIGIN   
   $!IF |Loop| == 1
      $!EXPORTSTART
        EXPORTREGION = ALLFRAMES
   $!ENDIF
   $!IF |Loop| != 1
      $!EXPORTNEXTFRAME
   $!ENDIF
$!ENDLOOP
$!EXPORTFINISH
