load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.debug.*.6.nc | head -n1) && echo ${a/\.debug.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + ".debug.*.6.nc")
; file_names = (/"mz_c_u_00.360x180.dt50.debug.0001-01-05T12:00:00Z.6.nc"/)

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

res_fv = res
res_fv@cnLevelSelectionMode = "ManualLevels"
res_fv@cnMinLevelValF = -0.3
res_fv@cnMaxLevelValF =  0.3
res_fv@cnLevelSpacingF = 0.01
res_fv@gsnLeftString = "fv"

res_fu = res
res_fu@cnLevelSelectionMode = "ManualLevels"
res_fu@cnMinLevelValF = -0.3
res_fu@cnMaxLevelValF =  0.3
res_fu@cnLevelSpacingF = 0.01
res_fu@gsnLeftString = "fu"

res_panel = True

do t = 0, dimsizes(file_names) - 1, 4
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  res_fv@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, f->fv(0,:,:), res_fv)
  printMinMax(f->fv, 0)
  res_fu@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot1 = gsn_csm_contour_map(wks, f->fu(0,:,:), res_fu)
  gsn_panel(wks, (/plot0, plot1/), (/2,1/), res_panel)
end do