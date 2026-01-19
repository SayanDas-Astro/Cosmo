import matplotlib.pyplot as plt
import numpy as np

try:
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    plt.plot(x, y)
    plt.savefig('test.png')
    print("SUCCESS")
except Exception as e:
    print(f"FAILED: {e}")
