import rasterio
import matplotlib.pyplot as plt
import numpy as np

# This code tests reading a .tif file, and generates a cool image with it!

with rasterio.open("mdt.tif") as src:
    print("Bounds:", src.bounds)
    print("CRS:", src.crs)
    print("NoData:", src.nodata)
    dem = src.read(1)

dem = np.where(dem == -9999, np.nan, dem)

valid = dem[~np.isnan(dem)]

print("Min:", np.min(valid))
print("Max:", np.max(valid))
print("Mean:", np.mean(valid))

plt.figure(figsize=(12,8))
plt.imshow(
    dem,
    vmin=np.percentile(valid, 2),
    vmax=np.percentile(valid, 98)
)
plt.colorbar(label="Altitude (m)")
plt.show()