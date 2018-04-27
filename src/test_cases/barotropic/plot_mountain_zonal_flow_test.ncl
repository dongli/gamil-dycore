if (.not. isdefined("file_prefix")) then
  system("echo '[Error]: Argument file_prefix should be set!'")
  exit
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")


wks = gsn_open_wks("pdf", file_prefix)

res_h = True
res_h@gsnDraw = False
res_h@gsnFrame = False
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 5050
res_h@cnMaxLevelValF = 5950
res_h@cnLevelSpacingF = 100
res_h@gsnPolar = "NH"

res_u = True
res_u@gsnDraw = False
res_u@gsnFrame = False
res_u@vcGlyphStyle = "CurlyVector"
res_u@gsnPolar = "NH"

res_panel = True
res_panel@gsnMaximize = True

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
  plot0 = gsn_csm_contour_map(wks, h, res_h)

  u = f->u(0,::5,::5)
  v = f->v(0,::5,::5)
  plot1 = gsn_csm_vector_map(wks, u, v, res_u)

  gsn_panel(wks, (/plot0,plot1/), (/2,1/), res_panel)
end do