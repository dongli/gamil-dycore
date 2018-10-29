load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.debug.*.3.nc | head -n1) && echo ${a/\.debug.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + ".debug.*.3.nc")

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

res_v_adv_lon = res
res_v_adv_lon@cnLevelSelectionMode = "ManualLevels"
res_v_adv_lon@cnMinLevelValF = -0.01
res_v_adv_lon@cnMaxLevelValF =  0.01
res_v_adv_lon@cnLevelSpacingF = 0.001
res_v_adv_lon@gsnLeftString = "v_adv_lon"

res_v_adv_lat = res
res_v_adv_lat@cnLevelSelectionMode = "ManualLevels"
res_v_adv_lat@cnMinLevelValF = -0.01
res_v_adv_lat@cnMaxLevelValF =  0.01
res_v_adv_lat@cnLevelSpacingF = 0.001
res_v_adv_lat@gsnLeftString = "v_adv_lat"

res_panel = True

do t = 0, dimsizes(file_names) - 1, 4
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  res_v_adv_lon@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, f->v_adv_lon(0,:,:), res_v_adv_lon)
  printMinMax(f->v_adv_lon, 0)
  res_v_adv_lat@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot1 = gsn_csm_contour_map(wks, f->v_adv_lat(0,:,:), res_v_adv_lat)
  gsn_panel(wks, (/plot0, plot1/), (/2,1/), res_panel)
end do