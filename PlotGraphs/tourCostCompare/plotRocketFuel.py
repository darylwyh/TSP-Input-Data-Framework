import matplotlib.pyplot as plt

# Concorde算法的Tour cost
concorde_costs = [2.59, 2.72, 2.81, 2.95, 3.32, 3.56]
concorde_sizes = [10, 20, 25, 30, 40, 50]

# 随机插入启发式算法的Tour cost
random_insertion_costs = [2.59, 2.72, 2.81, 2.82, 2.96, 3.33, 3.88]
random_insertion_sizes = [10, 20, 23, 25, 30, 40, 50]

# Held-Karp算法的Tour cost
held_karp_costs = [2.59, 2.72, 2.718, 2.735, 2.805, 2.808]
held_karp_sizes = [10, 20, 21, 22, 23, 24]

# Attention Model算法的Tour cost
attention_model_costs = [2.59, 2.72, 2.83, 2.95, 3.41, 3.61]
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
plt.title('RocketFuel Tour Cost Comparison of Different Algorithms')
plt.legend()
plt.grid(True)
plt.show()
