#!MC 1120
# Created by Tecplot 360 build 11.3.29.563
$!READDATASET  ' "../fort.100" "../fort.101" "../fort.102" "../fort.103" "../fort.104" "../fort.105" "../fort.106" "../fort.107" "../fort.108" "../fort.109" "../fort.110" "../fort.111" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "z" "p"'
$!GLOBALCONTOUR 1  VAR = 4
$!ISOSURFACELAYERS SHOW = YES
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 0

$!VarSet |first_line| = |NUMZONES|
$!VarSet |first_line| += 1
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'Extract Over Time'
  COMMAND = 'ExtractIsoSurfaceOverTime'
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