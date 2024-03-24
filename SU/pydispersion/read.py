import numpy as np
import matplotlib.pyplot as plt

def su(filename, endian='little'):
    """ Reads Seismic Unix files
    """
    import obspy
    if endian=='big':
        stream = obspy.read(filename, 
                   format='SU',
                   byteorder='>')
    else:
        stream = obspy.read(filename, 
                   format='SU',
                   byteorder='<')
    return stream

su_data = readers.su('/Users/zhiyuzhang/MyProjects/PyDispersion/data/viscoelastic_waves/Uz_file_single.su', endian='big')

print(len(su_data))                                                 # 50
data_matrix = np.vstack([trace.data for trace in su_data])          # (50, 10000)
data_matrix = np.column_stack([trace.data for trace in su_data])    # (10000, 50)
print(data_matrix.shape)

plt.imshow(data_matrix, aspect='auto', cmap='seismic', interpolation='bilinear', extent=[0, 50, 0, 10000])
