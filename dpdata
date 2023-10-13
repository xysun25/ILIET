import dpdata
import numpy as np

outcar_data = dpdata.LabeledSystem("OUTCAR")

coords = outcar_data["coords"]
forces = outcar_data["forces"]

# 获取三维数组的形状
original_shape = coords.shape
# 计算新的形状，将前两个维度合并为一个，保持第三个维度不变
new_shape = (original_shape[0] * original_shape[1], original_shape[2])
# 使用reshape方法将三维数组变为二维数组
two_coords = coords.reshape(new_shape)

# 获取三维数组的形状
original_shape1 = forces.shape
# 计算新的形状，将前两个维度合并为一个，保持第三个维度不变
new_shape1 = (original_shape1[0] * original_shape1[1], original_shape1[2])
# 使用reshape方法将三维数组变为二维数组
two_forces = forces.reshape(new_shape1)

result = np.concatenate((two_coords, two_forces), axis=1)

output_file = 'output.txt'
with open(output_file, 'w') as f:
    for row in result:
        f.write(' '.join(map(str, row)) + '\n')
