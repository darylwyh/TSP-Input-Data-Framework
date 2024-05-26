import matplotlib.pyplot as plt
# China数据集的Elapsed Time和Memory Usage
concorde_sizes = [10, 20, 25, 30, 40, 50]
concorde_times = [0.01, 0.02, 0.03, 0.04, 0.04, 0.04]
concorde_memory = [351.25, 352.04, 350.14, 351.20, 350.58, 350.70]

random_insertion_sizes = [10, 20, 25, 30, 40, 50]
random_insertion_times = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
random_insertion_memory = [351.80, 352.08, 351.61, 350.32, 350.75, 350.17]

held_karp_sizes = [20, 21, 22, 23, 24]
held_karp_times = [1.187, 3.296, 6.780, 15.569, 38.609]
held_karp_memory = [254, 493, 992, 2016, 3823]

attention_model_sizes = [10, 20, 25, 30, 40, 50]
attention_model_times = [0.01, 0.00, 0.01, 0.01, 0.01, 0.01]
attention_model_memory = [511.67, 514.59, 517.93, 520.84, 520.84, 520.84]

plt.figure(figsize=(10, 6))

# 绘制Elapsed Time
plt.subplot(2, 1, 1)
plt.plot(concorde_sizes, concorde_times, marker='o', label='Concorde')
plt.plot(random_insertion_sizes, random_insertion_times, marker='x', label='Random Insertion')
plt.plot(held_karp_sizes, held_karp_times, marker='s', label='Held-Karp')
plt.plot(attention_model_sizes, attention_model_times, marker='d', label='Attention Model')
plt.xlabel('Number of Cities (N)')
plt.ylabel('Elapsed Time (seconds)')
plt.title('China Elapsed Time Comparison of Different Algorithms')
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
plt.title('China Memory Usage Comparison of Different Algorithms')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
