"""Cut every sensor's image to the same square window (UTM 46N) around a validation point."""
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from pyproj import Transformer
CRS = 'EPSG:32646'
def window(lon, lat, half_m):
    x, y = Transformer.from_crs(4326, CRS, always_xy=True).transform(lon, lat)
    return x - half_m, y - half_m, x + half_m, y + half_m, (x, y)
def chip_rgb(path, win, px, bands=(1, 2, 3), resampling=Resampling.bilinear, scale=None, src_nodata=None):
    """Return float RGB (px x px) of `path` on the window; scale=None keeps 0-255 bytes -> 0-1."""
    x0, y0, x1, y1, _ = win
    tr = from_origin(x0, y1, (x1 - x0) / px, (y1 - y0) / px)
    out = np.full((3, px, px), np.nan, 'float32')
    with rasterio.open(path) as r:
        for i, b in enumerate(bands):
            src = r.read(b).astype('float32')
            if src_nodata is not None: src[src == src_nodata] = np.nan
            reproject(src, out[i], src_transform=r.transform, src_crs=r.crs, dst_transform=tr, dst_crs=CRS,
                      resampling=resampling, src_nodata=np.nan, dst_nodata=np.nan)
    rgb = np.moveaxis(out, 0, -1)
    return (rgb / 255.0 if scale is None else rgb / scale), tr
