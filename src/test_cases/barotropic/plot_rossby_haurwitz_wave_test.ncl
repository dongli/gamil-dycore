load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if
if (.not. isdefined("NH")) then
  NH = False
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

wks = gsn_open_wks("pdf", file_prefix)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@mpOutlineOn = False
res@mpFillOn = False
res@cnLineLabelsOn = False
if (NH) then
  res@gsnPolar = "NH"
end if

res_h = res
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 8000
res_h@cnMaxLevelValF = 11000
res_h@cnLevelSpacingF = 200
res_h@gsnLeftString = "Height (m)"

res_u = res
res_u@cnLevelSelectionMode = "ManualLevels"
res_u@cnMinLevelValF = -10
res_u@cnMaxLevelValF = 100
res_u@cnLevelSpacingF = 10
res_u@gsnLeftString = "U (m/s)"

res_v = res
res_v@cnLevelSelectionMode = "ManualLevels"
res_v@cnMinLevelValF = -70
res_v@cnMaxLevelValF = 70
res_v@cnLevelSpacingF = 10
res_v@gsnLeftString = "V (m/s)"


do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")

  time = f->time

  h = f->gh(0,:,:) / 9.80616
  h!0 = "lat"
  h!1 = "lon"
  h&lat = f->lat
  h&lon = f->lon
  h@long_name = "Height"
  h@units = "m"
  res_h@gsnRightString = cd_string(time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, h, res_h)

  u = f->u(0,:,:)
  res_u@gsnRightString = cd_string(time, "Day %j Hour %H")
  plot1 = gsn_csm_contour_map(wks, u, res_u)

  v = f->v(0,:,:)
  res_v@gsnRightString = cd_string(time, "Day %j Hour %H")
  plot2 = gsn_csm_contour_map(wks, v, res_v)

  gsn_panel(wks, (/plot0, plot1, plot2/), (/3, 1/), False)
end do
