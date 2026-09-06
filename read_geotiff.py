import rasterio
import matplotlib.pyplot as plt
import numpy as np

# This code tests reading a .tif file, and generates a cool image with it!

with rasterio.open("mdt.tif") as src:
    bounds = src.bounds
    print("Bounds:", src.bounds)
    print("CRS:", src.crs)
    print("NoData:", src.nodata)
    print("Dimensões (Linhas, Colunas):", src.shape)
    print("Resolução do Pixel (X, Y) em graus:", src.res)
    print("Matriz de Transformação Afim:\n", src.transform)

    print("----------------------------------------")
    print(f"CRS do Tile: {src.crs}")
    print(f"Limites (Esquerda, Fundo, Direita, Topo):")
    print(f"minx (Longitude mínima): {bounds.left}")
    print(f"miny (Latitude mínima): {bounds.bottom}")
    print(f"maxx (Longitude máxima): {bounds.right}")
    print(f"maxy (Latitude máxima): {bounds.top}")
    print("----------------------------------------")
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