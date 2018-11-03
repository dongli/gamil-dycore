load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if
if (.not. isdefined("NH")) then
  NH = False
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

wks = gsn_open_wks("pdf", file_prefix + ".hov")

gsn_define_colormap(wks, "BlWhRe")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@mpOutlineOn = False
res@mpFillOn = False
res@cnLineLabelsOn = False
res@cnFillOn = True
res@cnLinesOn = False
if (NH) then
  res@gsnPolar = "NH"
end if

res_h = res
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = -1
res_h@cnMaxLevelValF =  1
res_h@cnLevelSpacingF = 0.1
res_h@gsnLeftString = "Height (m)"
res_h@lbLabelBarOn = False

res_u = res
res_u@cnLevelSelectionMode = "ManualLevels"
res_u@cnMinLevelValF = -1
res_u@cnMaxLevelValF =  1
res_u@cnLevelSpacingF = 0.1
res_u@gsnLeftString = "U (m/s)"
res_u@lbLabelBarOn = False

res_v = res
res_v@cnLevelSelectionMode = "ManualLevels"
res_v@cnMinLevelValF = -1
res_v@cnMaxLevelValF =  1
res_v@cnLevelSpacingF = 0.1
res_v@gsnLeftString = "V (m/s)"
res_v@lbLabelBarOn = True

; days = 0
; do t = 0, dimsizes(file_names) - 1
;   system("echo '=> " + file_names(t) + "'")
;   f = addfile(file_names(t), "r")

;   time = f->time

;   h = f->gd(0,:,:) - 5.0e4
;   h = h / max(abs(h))
;   h!0 = "lat"
;   h!1 = "lon"
;   h&lat = f->lat
;   h&lon = f->lon
;   res_h@gsnRightString = "Days " + days
;   ; res_h@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot0 = gsn_csm_contour_map(wks, h, res_h)

;   u = f->u(0,:,:) / max(abs(f->u(0,:,:)))
;   u!0 = "lat"
;   u!1 = "lon"
;   u&lat = f->lat
;   u&lon = f->lon
;   res_u@gsnRightString = "Days " + days
;   ; res_u@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot1 = gsn_csm_contour_map(wks, u, res_u)

;   v = f->v(0,:,:) / max(abs(f->v(0,:,:)))
;   v!0 = "lat"
;   v!1 = "lon"
;   v&lat = f->lat
;   v&lon = f->lon
;   res_v@gsnRightString = "Days " + days
;   ; res_v@gsnRightString = cd_string(time, "Day %j Hour %H")
;   plot2 = gsn_csm_contour_map(wks, v, res_v)

;   gsn_panel(wks, (/plot0, plot1, plot2/), (/3, 1/), False)

;   days = days + 10
; end do

res_hov = True
res_hov@cnLevelSelectionMode = "ManualLevels"
res_hov@cnMinLevelValF = -1
res_hov@cnMaxLevelValF =  1
res_hov@cnLevelSpacingF = 0.1
res_hov@cnFillOn = True
res_hov@cnLinesOn = False
res_hov@cnLineLabelsOn = False
res_hov@cnFillMode = "RasterFill"
; res_h@gsnLeftString = "Height (m)"

f = addfile(file_names(0), "r")
v = new((/dimsizes(file_names), dimsizes(f->lon)/), typeof(f->v))
do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  v(t,:) = f->v(0,{44},:)
  v(t,:) = (/v(t,:) / max(abs(v(t,:)))/)
end do
v!0 = "time"
v!1 = "lon"
v&lon = f->lon
plot = gsn_csm_contour(wks, v, res_hov)
draw(plot)
frame(wks)