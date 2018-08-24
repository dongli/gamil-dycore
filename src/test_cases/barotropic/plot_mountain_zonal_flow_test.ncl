load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

wks = gsn_open_wks("pdf", file_prefix)

res_h = True
res_h@gsnDraw = False
res_h@gsnFrame = False
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 5050
res_h@cnMaxLevelValF = 5950
res_h@cnLevelSpacingF = 50
res_h@gsnPolar = "NH"
res_h@gsnLeftString = "Height (m)"

res_u = True
res_u@gsnDraw = False
res_u@gsnFrame = False
res_u@vcGlyphStyle = "CurlyVector"
res_u@gsnPolar = "NH"

res_panel = True
res_panel@gsnMaximize = True

res_h_zonal = True
res_h_zonal@trYMinF = 5025
res_h_zonal@trYMaxF = 5050

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
  res_h@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot0 = gsn_csm_contour_map(wks, h, res_h)
  draw(plot0)
  frame(wks)

  ;u = f->u(0,::5,::5)
  ;v = f->v(0,::5,::5)
  ;plot1 = gsn_csm_vector_map(wks, u, v, res_u)

  ;gsn_panel(wks, (/plot0,plot1/), (/2,1/), res_panel)

  ; res_h_zonal@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  ; plot = gsn_csm_xy(wks, f->lon, h(179,:), res_h_zonal)
end do
