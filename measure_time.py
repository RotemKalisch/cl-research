import time

def measure_time(func, *args):
    start = time.time()
    func(*args)
    end = time.time()
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(end - start)))
