load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

f = addfile(file_names(0), "r")
hs = f->ghs(0,:,:) / 9.80616

wks = gsn_open_wks("pdf", file_prefix)

res_h = True
res_h@gsnDraw = False
res_h@gsnFrame = False
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 5050
res_h@cnMaxLevelValF = 5950
res_h@cnLevelSpacingF = 50
res_h@gsnPolar = "NH"
res_h@mpOutlineOn = False
res_h@gsnLeftString = "Height (m)"

res_u = True
res_u@gsnDraw = False
res_u@gsnFrame = False
res_u@cnLevelSelectionMode = "ManualLevels"
res_u@cnMinLevelValF = -20
res_u@cnMaxLevelValF =  20
res_u@cnLevelSpacingF = 2
res_u@gsnPolar = "NH"
res_u@gsnLeftString = "U (m/s)"

res_v = True
res_v@gsnDraw = False
res_v@gsnFrame = False
res_v@cnLevelSelectionMode = "ManualLevels"
res_v@cnMinLevelValF = -20
res_v@cnMaxLevelValF =  20
res_v@cnLevelSpacingF = 2
res_v@gsnPolar = "NH"
res_v@gsnLeftString = "V (m/s)"

res_panel = True
res_panel@gsnMaximize = True
res_panel@gsnPanelRowSpec = True

do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  h = f->gh(0,:,:) / 9.80616
  h!0 = "lat"
  h!1 = "lon"
  h&lat = f->lat
  h&lon = f->lon
  res_h@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot0 = gsn_csm_contour_map(wks, h, res_h)

  res_u@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot1 = gsn_csm_contour_map(wks, f->u(0,:,:), res_u)

  res_v@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot2 = gsn_csm_contour_map(wks, f->v(0,:,:), res_v)

  gsn_panel(wks, (/plot0,plot1,plot2/), (/1,2/), res_panel)
end do

; res_h@cnFillOn = True
; res_h@cnFillPalette = "MPL_viridis"
; res_h@lbBoxEndCapStyle = "TriangleBothEnds"
; res_h@pmLabelBarHeightF = 0.05
; 
; res_uv_vector = True
; res_uv_vector@gsnDraw = False
; res_uv_vector@gsnFrame = False
; res_uv_vector@vcGlyphStyle = "LineArrow"
; res_uv_vector@vcRefMagnitudeF = 10
; res_uv_vector@vcRefLengthF = 0.01
; res_uv_vector@gsnLeftString = ""
; res_uv_vector@gsnRightString = ""
; 
; uv_stride = 5
; 
; do t = 0, dimsizes(file_names) - 1
;   system("echo '=> " + file_names(t) + "'")
;   f = addfile(file_names(t), "r")
;   d = f->gh(0,:,:) / 9.80616
;   h = hs + d
;   h!0 = "lat"
;   h!1 = "lon"
;   h&lat = f->lat
;   h&lon = f->lon
;   res_h@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
;   plot0 = gsn_csm_contour_map(wks, h, res_h)
; 
;   u := f->u(0,::uv_stride,::uv_stride)
;   v := f->v(0,::uv_stride,::uv_stride)
;   plot1 = gsn_csm_vector(wks, u, v, res_uv_vector)
; 
;   overlay(plot0, plot1)
; 
;   draw(plot0)
;   frame(wks)
; end do

; res_hs = True
; res_hs@mpOutlineOn = False
; res_hs@gsnLeftString = "Mountain (m)"

; f = addfile(file_names(0), "r")
; hs = f->ghs(0,:,:) / 9.80616
; hs!0 = "lat"
; hs!1 = "lon"
; hs&lat = f->lat
; hs&lon = f->lon
; plot = gsn_csm_contour_map(wks, hs, res_hs)

; res_hs := True
; res_hs@gsnLeftString = "Mountain (m)"
; res_hs@xyLineThicknessF = 5.0
; res_hs@trYMinF = -10

; plot = gsn_csm_xy(wks, (/f->lon/), hs(118,:), res_hs)
