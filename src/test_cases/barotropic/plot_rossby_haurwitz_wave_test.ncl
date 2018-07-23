load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  system("echo '[Error]: Argument file_prefix should be set!'")
  exit
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")


wks = gsn_open_wks("pdf", file_prefix)

res_h = True
res_h@gsnDraw = False
res_h@gsnFrame = False
res_h@mpOutlineOn = False
res_h@mpFillOn = False
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 8000
res_h@cnMaxLevelValF = 11000
res_h@cnLevelSpacingF = 200
res_h@cnLineLabelsOn = False
;res_h@gsnPolar = "NH"
res_h@gsnLeftString = "Height (m)"

do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  h = f->gd(0,:,:) / 9.80616
  h!0 = "lat"
  h!1 = "lon"
  h&lat = f->lat
  h&lon = f->lon
  h@long_name = "Height"
  h@units = "m"
  res_h@gsnRightString = cd_string(f->time, "Day %j Hour %H")
  plot0 = gsn_csm_contour_map(wks, h, res_h)
  draw(plot0)
  frame(wks)
end do
