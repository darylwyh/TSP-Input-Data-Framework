import matplotlib.pyplot as plt

# Berlin数据集的Elapsed Time和Memory Usage
concorde_sizes = [10, 20, 25, 30, 40, 50]
concorde_times = [0.01, 0.02, 0.03, 0.05, 0.05, 0.09]
concorde_memory = [350.16, 351.57, 350.77, 351.06, 351.57, 352.07]

random_insertion_sizes = [10, 20, 25, 30, 40, 50]
random_insertion_times = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
random_insertion_memory = [351.54, 350.09, 350.95, 351.23, 351.13, 350.10]

held_karp_sizes = [10, 20, 21, 22, 23, 24]
held_karp_times = [0.049, 1.171, 2.712, 6.036, 15.569, 37.628]
held_karp_memory = [2, 255, 486, 990, 2011, 3825]

attention_model_sizes = [10, 20, 25, 30, 40, 50]
attention_model_times = [0.00, 0.00, 0.00, 0.00, 0.01, 0.00]
attention_model_memory = [485.54, 489.77, 492.93, 497.34, 501.47, 508.45]

plt.figure(figsize=(10, 6))

# 绘制Elapsed Time
plt.subplot(2, 1, 1)
plt.plot(concorde_sizes, concorde_times, marker='o', label='Concorde')
plt.plot(random_insertion_sizes, random_insertion_times, marker='x', label='Random Insertion')
plt.plot(held_karp_sizes, held_karp_times, marker='s', label='Held-Karp')
plt.plot(attention_model_sizes, attention_model_times, marker='d', label='Attention Model')
plt.xlabel('Number of Cities (N)')
plt.ylabel('Elapsed Time (seconds)')
plt.title('Berlin Elapsed Time Comparison of Different Algorithms')
plt.legend()
plt.grid(True)

# 绘制Memory Usage
plt.subplot(2, 1, 2)
plt.plot(concorde_sizes, concorde_memory, marker='o', label='Concorde')
plt.plot(random_insertion_sizes, random_insertion_memory, marker='x', label='Random Insertion')
plt.plot(held_karp_sizes, held_karp_memory, marker='s', label='Held-Karp')
plt.plot(attention_model_sizes, attention_model_memory, marker='d', label='Attention Model')
plt.xlabel('Number of Cities (N)')
plt.ylabel('Memory Usage (MB)')
plt.title('Berlin Memory Usage Comparison of Different Algorithms')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
