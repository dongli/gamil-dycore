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

res_mass_div_lon = res
res_mass_div_lon@cnLevelSelectionMode = "ManualLevels"
res_mass_div_lon@cnMinLevelValF = -0.3
res_mass_div_lon@cnMaxLevelValF =  0.3
res_mass_div_lon@cnLevelSpacingF = 0.01
res_mass_div_lon@gsnLeftString = "mass_div_lon"

res_mass_div_lat = res
res_mass_div_lat@cnLevelSelectionMode = "ManualLevels"
res_mass_div_lat@cnMinLevelValF = -0.3
res_mass_div_lat@cnMaxLevelValF =  0.3
res_mass_div_lat@cnLevelSpacingF = 0.01
res_mass_div_lat@gsnLeftString = "mass_div_lat"

res_panel = True

do t = 0, dimsizes(file_names) - 1, 4
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  res_mass_div_lon@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, f->mass_div_lon(0,:,:), res_mass_div_lon)
  printMinMax(f->mass_div_lon, 0)
  res_mass_div_lat@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot1 = gsn_csm_contour_map(wks, f->mass_div_lat(0,:,:), res_mass_div_lat)
  gsn_panel(wks, (/plot0, plot1/), (/2,1/), res_panel)
end do