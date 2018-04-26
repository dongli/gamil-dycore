if (.not. isdefined("file_prefix")) then
  system("echo '[Error]: Argument file_prefix should be set!'")
  exit
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")


wks = gsn_open_wks("pdf", file_prefix)

res = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5050
res@cnMaxLevelValF = 5950
res@cnLevelSpacingF = 100
res@gsnPolar = "NH"

f = addfile(file_names(0), "r")
hs = f->ghs(0,:,:) / 9.80616
do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  d = f->gd(0,:,:) / 9.80616
  h = hs + d
  h!0 = "lat"
  h!1 = "lon"
  h&lat = f->lat
  h&lon = f->lon
  h@long_name = "Height"
  h@units = "m"
  plot = gsn_csm_contour_map(wks, h, res)
end do