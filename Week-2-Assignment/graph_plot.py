import matplotlib.pyplot as plt

fig, ax = plt.subplots()

probabilities = ['Ideal 0', '0.001', '0.005', '0.01', '0.05', '0.1']
counts = [0, 0, 0, 2, 5, 7]
percantage = ['0%', '0%', '0%', '20%', '50%', '70%']
bar_labels = ['Ideal 0', '0.001', '0.005', '0.01', '0.05', '0.1']
bar_colors = ['tab:purple' ,'tab:red', 'tab:blue', 'tab:red', 'tab:orange', 'tab:green']

ax.bar(probabilities, counts, label=bar_labels, color=bar_colors)

ax.set_ylabel('Number of logical errors in 10 Runs')
ax.set_xlabel('Probability of physical errors')
ax.set_title('Physical Error Probability vs Logical Bit-Flip Errors')
ax.legend(title='Physical Error Probability')

ax.bar_label(ax.containers[0], percantage ,padding=3, fontsize=15, fmt='%d')

plt.show()