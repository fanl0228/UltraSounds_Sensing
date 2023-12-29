import numpy as np
import matplotlib.pyplot as plt

def one_dimensional_kalman_filter(measurements, initial_state_mean, initial_state_covariance, process_variance, measurement_variance):
    state_mean = initial_state_mean
    state_covariance = initial_state_covariance

    filtered_states = []

    for measurement in measurements:
        # Prediction step
        state_mean_prior = state_mean
        state_covariance_prior = state_covariance + process_variance

        # Update step
        kalman_gain = state_covariance_prior / (state_covariance_prior + measurement_variance)
        state_mean = state_mean_prior + kalman_gain * (measurement - state_mean_prior)
        state_covariance = (1 - kalman_gain) * state_covariance_prior

        filtered_states.append(state_mean)

    return np.array(filtered_states)

if __name__ == '__main__':
    # 生成示例信号
    np.random.seed(42)
    true_signal = np.linspace(0, 10, 100) + np.random.normal(0, 1, 100)

    # 添加噪声作为测量值
    measurements = true_signal + np.random.normal(0, 0.5, 100)

    # 初始化滤波器参数
    initial_state_mean = 0
    initial_state_covariance = 1
    process_variance = 0.1
    measurement_variance = 0.5

    # 使用卡尔曼滤波
    filtered_states = one_dimensional_kalman_filter(measurements, initial_state_mean, initial_state_covariance, process_variance, measurement_variance)

    # 绘制结果
    plt.plot(true_signal, label='True Signal')
    plt.scatter(range(len(measurements)), measurements, color='red', marker='.', label='Measurements')
    plt.plot(filtered_states, label='Filtered Signal', color='green')
    plt.legend()
    plt.title('One-Dimensional Kalman Filter for Signal Estimation')
    plt.show()

    pass
