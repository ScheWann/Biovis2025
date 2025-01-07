from flask import Flask, request, send_file, jsonify
from PIL import Image
import os
import hashlib

app = Flask(__name__)

# 全局配置
IMAGE_PATH = '../Data/spatial/tissue_hires_images.png'  # 大图路径
CACHE_DIR = '../cache/'  # 缓存目录

# 确保缓存目录存在
os.makedirs(CACHE_DIR, exist_ok=True)

@app.route('/', methods=['GET'])
def get_helloword():
    return 'Hello World!'

@app.route('/get_tile', methods=['GET'])
def get_tile():
    try:
        # 获取请求参数
        x = int(request.args.get('x'))  # 瓦片 x 坐标
        y = int(request.args.get('y'))  # 瓦片 y 坐标
        z = int(request.args.get('z'))  # 缩放级别
        tile_size = int(request.args.get('tile_size', 256))  # 瓦片大小（默认 256x256）

        # 根据请求参数生成缓存键
        cache_key = f"{x}_{y}_{z}_{tile_size}"
        cache_path = os.path.join(CACHE_DIR, f"{cache_key}.png")

        # 如果缓存存在，直接返回
        if os.path.exists(cache_path):
            return send_file(cache_path, mimetype='image/png')

        # 加载原始图片
        image = Image.open(IMAGE_PATH)

        # 计算缩放比例和瓦片范围
        scale = 2 ** z  # 假设原图对应 z=0
        tile_x_min = x * tile_size // scale
        tile_y_min = y * tile_size // scale
        tile_x_max = (x + 1) * tile_size // scale
        tile_y_max = (y + 1) * tile_size // scale

        # 裁剪图片
        cropped_image = image.crop((tile_x_min, tile_y_min, tile_x_max, tile_y_max))
        cropped_image = cropped_image.resize((tile_size, tile_size), Image.Resampling.LANCZOS)

        # 保存到缓存
        cropped_image.save(cache_path, format='PNG')

        # 返回瓦片
        return send_file(cache_path, mimetype='image/png')

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)
