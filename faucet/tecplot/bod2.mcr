#!MC 1120
# Created by Tecplot 360 build 11.3.29.563

$!Varset |blocks| = 12
$!Varset |zero| = 100

$!NEWLAYOUT 
$!Varset |current| = |zero|

$!Varset |current| -= 1
$!Loop |blocks|
$!Varset |current| += 1

$!READDATASET  ' "fort.|current|.plt" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "p"'
$!Endloop

$!GLOBALCONTOUR 1  VAR = 4
$!ISOSURFACELAYERS SHOW = YES
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 0

$!VarSet |first_line| = |NUMZONES|
$!VarSet |first_line| += 1

# $!EXTENDEDCOMMAND 
#   COMMANDPROCESSORID = 'Extract Over Time'
#   COMMAND = 'ExtractIsoSurfaceOverTime'

$!CREATEISOZONES 

$!WRITEDATASET  "./bod.plt"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  ASSOCIATELAYOUTWITHDATAFILE = NO
  ZONELIST =  [|first_line|-|NUMZONES|]
  VARPOSITIONLIST =  [1-3]
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9