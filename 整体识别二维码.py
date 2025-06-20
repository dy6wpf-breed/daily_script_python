import cv2
from pyzbar.pyzbar import decode

def scan_qr_codes(image_path):
    # 读取图片
    image = cv2.imread(image_path)
    
    if image is None:
        print(f"无法读取图片: {image_path}")
        return 0
        
    # 解码二维码
    decoded_objects = decode(image)
    
    # 打印每个二维码的内容
    for obj in decoded_objects:
        print(f"二维码内容: {obj.data.decode('utf-8')}")
    
    return len(decoded_objects)

# 使用正斜杠的路径
num_qrcodes = scan_qr_codes('E:/desk/qr_code_image.jpg')

print(f"识别到 {num_qrcodes} 个二维码")