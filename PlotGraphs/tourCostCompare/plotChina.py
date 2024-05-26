import matplotlib.pyplot as plt

# China数据集的Tour cost
concorde_costs = [3.17, 4.16, 4.46, 4.58, 5.00, 5.37]
concorde_sizes = [10, 20, 25, 30, 40, 50]

random_insertion_costs = [3.54, 4.99, 5.60, 5.53, 6.06, 6.35]
random_insertion_sizes = [10, 20, 25, 30, 40, 50]

held_karp_costs = [4.16, 4.16, 4.30, 4.31, 4.27]
held_karp_sizes = [20, 21, 22, 23, 24]

attention_model_costs = [3.17, 4.17, 4.46, 4.58, 5.18, 5.61]
attention_model_sizes = [10, 20, 25, 30, 40, 50]

plt.figure(figsize=(10, 6))

# 绘制Concorde算法的Tour cost
plt.plot(concorde_sizes, concorde_costs, marker='o', label='Concorde')

# 绘制随机插入启发式算法的Tour cost
plt.plot(random_insertion_sizes, random_insertion_costs, marker='x', label='Random Insertion')

# 绘制Held-Karp算法的Tour cost
plt.plot(held_karp_sizes, held_karp_costs, marker='s', label='Held-Karp')

# 绘制Attention Model算法的Tour cost
plt.plot(attention_model_sizes, attention_model_costs, marker='d', label='Attention Model')

plt.xlabel('Number of Cities (N)')
plt.ylabel('Tour Cost')
plt.title('China Tour Cost Comparison of Different Algorithms')
plt.legend()
plt.grid(True)
plt.show()
