import matplotlib.pyplot as plt
# RocketFuel数据集的Elapsed Time和Memory Usage
concorde_sizes = [10, 20, 25, 30, 40, 50]
concorde_times = [0.03, 0.03, 0.03, 0.04, 0.06, 0.07]
concorde_memory = [351.53, 351.22, 350.77, 351.24, 351.48, 350.58]

random_insertion_sizes = [10, 20, 23, 25, 30, 40, 50]
random_insertion_times = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
random_insertion_memory = [351.29, 352.01, 350.07, 350.52, 350.86, 350.52, 350.50]

held_karp_sizes = [10, 20, 21, 22, 23, 24]
held_karp_times = [0.04, 1.07, 2.58, 6.24, 15.14, 39.82]
held_karp_memory = [2, 255, 490, 985, 2008, 3822]

attention_model_sizes = [10, 20, 25, 30, 40, 50]
attention_model_times = [0.01, 0.00, 0.00, 0.00, 0.01, 0.00]
attention_model_memory = [468.85, 468.85, 471.70, 475.56, 478.14, 482.62]

plt.figure(figsize=(10, 6))

# 绘制Elapsed Time
plt.subplot(2, 1, 1)
plt.plot(concorde_sizes, concorde_times, marker='o', label='Concorde')
plt.plot(random_insertion_sizes, random_insertion_times, marker='x', label='Random Insertion')
plt.plot(held_karp_sizes, held_karp_times, marker='s', label='Held-Karp')
plt.plot(attention_model_sizes, attention_model_times, marker='d', label='Attention Model')
plt.xlabel('Number of Cities (N)')
plt.ylabel('Elapsed Time (seconds)')
plt.title('RocketFuel Elapsed Time Comparison of Different Algorithms')
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
plt.title('RocketFuel Memory Usage Comparison of Different Algorithms')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
