import matplotlib.pyplot as plt

# Berlin数据集的Tour cost
concorde_costs = [2.96, 3.91, 4.08, 4.47, 4.92, 5.43]
concorde_sizes = [10, 20, 25, 30, 40, 50]

random_insertion_costs = [2.96, 4.13, 4.37, 4.79, 5.19, 5.83]
random_insertion_sizes = [10, 20, 25, 30, 40, 50]

held_karp_costs = [2.96, 3.91, 3.94, 3.94, 3.99, 4.03]
held_karp_sizes = [10, 20, 21, 22, 23, 24]

attention_model_costs = [2.96, 3.91, 4.08, 4.49, 5.28, 5.94]
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
plt.title('Berlin Tour Cost Comparison of Different Algorithms')
plt.legend()
plt.grid(True)
plt.show()
