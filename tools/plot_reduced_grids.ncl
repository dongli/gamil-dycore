dlon = 2.0
dlat = 2.0

rf = new(90, integer)
rf(:) = 1
rf(0:2) = (/4, 2, 2/)

wks = gsn_open_wks("pdf", "reduce_grids")

res_base = True
res_base@gsnDraw = False
res_base@gsnFrame = False
res_base@mpProjection = "Satellite"
res_base@mpCenterLonF = 105
res_base@mpCenterLatF = 45
res_base@mpPerimOn = False
res_base@mpGridAndLimbOn = True
res_base@mpGridLonSpacingF = dlon
res_base@mpGridLatSpacingF = dlat
res_base@mpGridLineColor = "gray"
res_base@mpGridLineThicknessF = 0.00005

res_reduce = True
res_reduce@gsFillColor = "green"
res_reduce@gsLineColor = "red"
res_reduce@gsLineThicknessF = 0.4

plot = gsn_csm_map(wks, res_base)

lon = new(5, double)
lat = new(5, double)
do j = 1, dimsizes(rf) - 2
  if (rf(j) .eq. 1) then
    continue
  end if
  do i = 0, toint(360 / dlon / rf(j)) - 1
    lon(0) = i * dlon * rf(j)
    lat(0) = 90.0 + dlat * 0.5 - j * dlat
    lon(1) = lon(0)
    lat(1) = lat(0) - dlat
    lon(2) = lon(1) + dlon * rf(j)
    lat(2) = lat(1)
    lon(3) = lon(2)
    lat(3) = lat(0)
    lon(4) = lon(0)
    lat(4) = lat(0)
    ; plot@$unique_string("reduce")$ = gsn_add_polygon(wks, plot, lon, lat, res_reduce)
    plot@$unique_string("reduce")$ = gsn_add_polyline(wks, plot, lon, lat, res_reduce)
  end do
end do

draw(plot)
frame(wks)
