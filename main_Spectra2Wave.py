# 2025年1月8日 17:17
'''
python环境3.11.5

将波浪谱转化为波高的时间序列

2025-01-08: 
波浪谱为JONSWAP谱
References:
[1] Mazzaretto O M, Menéndez M, Lobeto H. A global evaluation of the JONSWAP spectra suitability on coastal areas[J]. Ocean Engineering, 2022, 266: 112756.
[2] https://ww2.mathworks.cn/matlabcentral/answers/58194-how-to-generate-a-time-signal-from-spectrum * 主要参考matlab代码
[2] https://ww2.mathworks.cn/matlabcentral/answers/58194-how-to-generate-a-time-signal-from-spectrum * 主要参考matlab代码
[3] Goda Y. Random seas and design of maritime structures[M]. World scientific, 2010.

'''

import numpy as np
import matplotlib.pyplot  as plt

class class_wave:
    def __init__(self, wave_height, wave_period, gamma=3.3):    # 初始化波浪对象
        self.Hs = wave_height
        self.Tp = wave_period
        self.gamma = gamma

    def jonswap(self, frequence: np.array):   # 计算jonswap谱
        # 波浪谱参数
        Hs = self.Hs
        Tp = self.Tp
        gamma = self.gamma
        beta = 0.0624/(0.23 + 0.0336*gamma - 0.185/(1.9+gamma) ) * \
            (1.094 - 0.01915*np.log(gamma))
        # 计算sigma
        condition = frequence <= 1/Tp
        sigma = np.where(condition, 0.07, 0.09)

        # 计算波浪谱
        Spec = beta * Hs**2 * Tp**(-4) * frequence**(-5) * \
            np.exp(-1.25 * (Tp*frequence)**(-4)) * \
            gamma**(np.exp(-(Tp*frequence-1)**2 / (2 * sigma**2)))
        return Spec
    

    def spectra2wave(self, Spec: np.array, frequence: np.array, time_s: np.array):
        frequence_step = frequence[1]-frequence[0]
        # 波高的幅值
        amplititude = np.sqrt(2 * Spec * frequence_step)  # 一维数组
        # 波浪分量的数量
        wave_num = len(frequence)
        # 相位
        phase = np.random.rand(wave_num) * 2 * np.pi    # 一维数组
        phase = np.random.rand(wave_num) * 2 * np.pi    # 一维数组
        # 圆频率
        omega = frequence * 2 * np.pi   # 一维数组   # 一维数组

        # 波浪分量
        wave_componant = [[] for _ in range(wave_num)]  # 二维数组，每一行是一个波浪分量，每一列是一个时刻
        wave_componant = [[] for _ in range(wave_num)]  # 二维数组，每一行是一个波浪分量，每一列是一个时刻
        for i in range(wave_num):
            wave_componant[i] = amplititude[i] * np.cos(omega[i]*time_s + phase[i])

        # 波浪=分量之和
        series = np.sum(wave_componant, axis=0)
        
        return series
    
    def wave2spectra(self, wave_series, smooth_sample = 20):
        # 计算幅值
        amplititude = np.fft.fft(wave_series)
        # 计算幅值对应频率（自动划分）
        frequence = abs(np.fft.fftfreq(len(wave_series)))
        # 计算频率步长
        frequence_step = frequence[1] - frequence[0]
        # 计算功率谱
        Spec = amplititude.real**2 / 2 * frequence_step *2*np.pi
        
        # 谱平滑
        window_size = len(amplititude)//smooth_sample
        Spec_smooth = []
        for i in range(len(Spec)):
            # 计算滑动窗口的起始和结束索引
            start = max(0, i - window_size // 2)
            end = min(len(Spec), i + window_size // 2 + 1)

            # 计算滑动窗口内数据的平均值
            window_data = Spec[start:end]
            window_mean = sum(window_data) / len(window_data)
            Spec_smooth.append(window_mean)
            
        return Spec_smooth, frequence
        
    def draw_spectra(self, frequence: np.array, Spec: np.array):
        plt.plot(frequence, Spec)
        plt.xlabel('Frequence (Hz)')
        plt.ylabel('PSD (m**2*s)')
        plt.grid(True, which='major',axis='both')
        plt.tick_params(which='both', direction='in')  # 设置主刻度和次刻度的方向为内向
        plt.xlim([0, max(frequence)])
        plt.show()

    def draw_series(self, time_s: np.array, series: np.array):
        plt.plot(time_s, series)
        plt.xlabel('Time (s)')
        plt.ylabel('Wave height (m)')
        plt.grid(True, which='major',axis='both')
        plt.tick_params(which='both', direction='in')  # 设置主刻度和次刻度的方向为内向
        plt.xlim([0, max(time_s)])
        plt.show()



if __name__ == "__main__":
    # 波浪谱参数
    Hs = 2
    Tp = 5
    gamma = 3.3
    # 波浪谱频率范围
    frequence_max = 1/Tp * 3
    frequence_step = 0.0001
    f = np.arange(frequence_step, frequence_max, frequence_step)
    # 时域时间
    time_length = 3600
    time_step = 1
    t_series = np.arange(0, time_length, time_step)

    # 创建波浪对象
    wave = class_wave(Hs, Tp, gamma)
    # 创建jonswap谱
    wave_S = wave.jonswap(f)
    # 展示波浪谱
    wave.draw_spectra(f, wave_S)
    # 波浪时间序列
    wave_series = wave.spectra2wave(wave_S, f, t_series)
    # 绘制时间序列
    wave.draw_series(t_series, wave_series)
    
    # 波浪序列转化为波浪谱
    Spectrum, freqs = wave.wave2spectra(wave_series, smooth_sample=150)
    # 展示波浪谱
    wave.draw_spectra(freqs, Spectrum)
    # 比较理论jonswap谱和转化谱
    plt.plot(f, wave_S, 'b--')
    plt.plot(freqs, Spectrum, 'r-')
    plt.xlabel('Frequence (Hz)')
    plt.ylabel('Wave spectra density (m**2*s)')
    plt.legend(['Theoretical', 'Empirical'])
    plt.grid(True, which='major',axis='both')
    plt.tick_params(which='both', direction='in')  # 设置主刻度和次刻度的方向为内向
    plt.show()
    





    print('FIN!!')

