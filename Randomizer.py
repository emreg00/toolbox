
import random

def main():
    number_of_samples = 50000    
    items = range(1,10000)
    szie = 100
    rparser = Randomizer.Randomizer(items, size, number_of_samples)
    while True:
	try:
	    rparser.get_next_sample()
	except StopIteration:
	    break
    return

def get_random_sub_sample(population, sample_size, seed_value=None):
    if seed_value is not None:
	random.seed(seed_value)
    return random.sample(population, sample_size)

class Randomizer(object):

    def __init__(self, population, sample_size, number_of_samples=50000):
	#self.population = population
	#self.sample_size = sample_size
	#self.number_of_samples = number_of_samples
	self.generator = self.get_generator(population, sample_size, number_of_samples)

    def get_generator(self, population, sample_size, number_of_samples):
	for i in xrange(0, number_of_samples):
	    sub_sample = get_random_sub_sample(population, sample_size)
	    yield sub_sample

    def get_next_sample(self):
	self.current_sample = self.generator.next()
	return self.current_sample

    def get_current_sample(self):
	return self.current_sample

    def fetch_next(self):
	return self.get_next_sample()

    def get_current(self, state_idx=None):
	return self.current_sample


