if (.not. isdefined("file_path")) then
  system("echo '[Error]: Argument file_path should be set!'")
  exit
end if

f = addfile(file_path, "r")

dlon = 1.0
dlat = 1.0

rf = f->rf(0,:)

wks = gsn_open_wks("pdf", "reduce_grids")

res_base = True
res_base@gsnDraw = False
res_base@gsnFrame = False
res_base@mpProjection = "Satellite"
res_base@mpCenterLonF = 0
res_base@mpCenterLatF = 45
res_base@mpPerimOn = False
res_base@mpGridAndLimbOn = True
res_base@mpGridLonSpacingF = 1
res_base@mpGridLatSpacingF = 1
res_base@mpGridLineThicknessF = 0.0001

res_reduce = True
res_reduce@gsFillColor = "green"
res_reduce@gsLineColor = "red"

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