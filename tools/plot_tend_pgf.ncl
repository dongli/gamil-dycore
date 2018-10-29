load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.debug.*.6.nc | head -n1) && echo ${a/\.debug.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + ".debug.*.6.nc")

wks = gsn_open_wks("pdf", file_prefix)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@mpOutlineOn = False
res@mpFillOn = False
res@cnConstFEnableFill = True
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnFillOn = True
res@cnFillMode = "RasterFill"
res@cnFillPalette = "amwg_blueyellowred"

res_u_pgf = res
res_u_pgf@cnLevelSelectionMode = "ManualLevels"
res_u_pgf@cnMinLevelValF = -0.5
res_u_pgf@cnMaxLevelValF =  0.5
res_u_pgf@cnLevelSpacingF = 0.01
res_u_pgf@gsnLeftString = "u_pgf"

res_v_pgf = res
res_v_pgf@cnLevelSelectionMode = "ManualLevels"
res_v_pgf@cnMinLevelValF = -0.5
res_v_pgf@cnMaxLevelValF =  0.5
res_v_pgf@cnLevelSpacingF = 0.01
res_v_pgf@gsnLeftString = "v_pgf"

res_panel = True

do t = 0, dimsizes(file_names) - 1, 4
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  res_u_pgf@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, f->u_pgf(0,:,:), res_u_pgf)
  printMinMax(f->u_pgf, 0)
  res_v_pgf@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot1 = gsn_csm_contour_map(wks, f->v_pgf(0,:,:), res_v_pgf)
  gsn_panel(wks, (/plot0, plot1/), (/2,1/), res_panel)
end do