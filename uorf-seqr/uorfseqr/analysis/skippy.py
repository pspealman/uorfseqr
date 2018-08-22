#Copyright MIT License, 2017 Armaghan Naik, Pieter Spealman

import numpy
# numpy.random.seed(4)

class OverInsertionError(Exception):
	pass

class Skipnode(object):
	def __init__(self, height, value, nexts):
		self.height = height
		self.next = nexts
		self.value = value
		self.markers = [set([]) for x in range(len(nexts))]
		self.endpoint_markers = set([])
		self.starting_markers = set([])
		self.ending_markers = set([])


class IntervalSkiplist(object):
	def __init__(self, maxheight):
		self.maxheight = maxheight
		self.tail = Skipnode(maxheight, numpy.inf, [None] * maxheight)
		self.head = Skipnode(maxheight, -numpy.inf, [self.tail] * maxheight)
		self.lookup = {}
		self.distribution = (1-(0.25)**numpy.arange(1,maxheight+1))
		self.distribution /= numpy.max(self.distribution)
		# print 'dist', self.distribution

	def horizon(self, value):
		atlevel = self.maxheight - 1
		place = self.head
		frontier = [None] * self.maxheight
		while (atlevel >= 0):
			if place.next[atlevel].value < numpy.inf and place.next[atlevel].value < value:
				place = place.next[atlevel]
			else:
				frontier[atlevel] = place
				atlevel -= 1
		# print 'frontier is', [x.value for x in frontier if x is not None]
		return frontier

	def insert_node(self, value):
		frontier = self.horizon(value)
		closest_node = frontier[0].next[0] 
		if closest_node.value == value:
			# already in skiplist
			return closest_node
		r = numpy.random.rand(1)[0]
		# print r, self.distribution
		v = numpy.where(r <= self.distribution )[0][-1]+1
		nexts = [None] * self.maxheight
		me = Skipnode(v, value, nexts)
		for h in range(0, v):
			me.next[h] = frontier[h].next[h]
			frontier[h].next[h] = me
		self.adjust_markers_after_insert(me, frontier)
		# print '------INSERTED', value, '@', v
		# self.printstuff()
		return me

	def adjust_markers_after_insert(self, node, frontier):
		promoted = []
		# print 'adjusting marker', node.value
		# Add the markers leading out of node at the highest possible level
		for alevel in range(node.height-2):
			newpromoted = []
			for marker in frontier[alevel].markers[alevel].copy():
				# print 'considering marker', marker, 'at', alevel
				starti, endi = self.lookup[marker]
				# if the interval of this marker contains (node.value, node.next[alevel+1].value)
				# print 'at node', node.value, 'at level', alevel, node.next[alevel+1].value
				if node.next[alevel + 1].value <= endi:
					# promote this marker
					self.remove_marker_on_path(marker, node.next[alevel], node.next[alevel + 1], alevel)
					newpromoted.append(marker)
				else:
					# place the marker on the level atlevel edge out of node
					node.markers[alevel].add(marker)
			for marker in promoted:
				starti, endi = self.lookup[marker]
				# if the interval of marker does *not* contain (node.value, node.next[alevel+1].value)
				# then we don't need to promote higher (reverse cases from paper)
				if node.next[alevel + 1].value < endi:
					# continue to promote this marker
					self.remove_marker_on_path(marker, node.next[alevel], node.next[alevel + 1], alevel)
					newpromoted.append(marker)
				else:
					# drop this marker off here
					node.markers[alevel].add(marker)
			promoted = newpromoted
		# finish off
		alevel = node.height-1
		node.markers[alevel] |= set(promoted) | frontier[alevel].markers[alevel]
		# now look at the markers leading into node and push them into levels at most this node's height
		promoted = []
		# Add the markers leading out of node at the highest possible level
		for alevel in range(node.height-2):
			newpromoted = []
			for marker in frontier[alevel].markers[alevel].copy():
				starti, endi = self.lookup[marker]
				# if the interval of this marker contains (node.value, node.next[alevel+1].value)
				if starti < frontier[alevel + 1].value:
					# promote this marker
					self.remove_marker_on_path(marker, frontier[alevel + 1], node, alevel)
					newpromoted.append(marker)
			for marker in promoted:
				starti, endi = self.lookup[marker]
				# if the interval of marker does *not* contain (node.value, node.next[alevel+1].value)
				# then we don't need to promote higher (reverse cases from paper)
				if starti < frontier[alevel + 1].value:
					# continue to promote this marker
					self.remove_marker_on_path(marker, frontier[alevel + 1], node, alevel)
					newpromoted.append(marker)
				else:
					frontier[alevel].markers[alevel].add(marker)
			promoted = newpromoted
		# finish off
		alevel = node.height-1
		frontier[alevel].markers[alevel] |= set(promoted) | frontier[alevel].markers[alevel]
		return

	def remove_marker_on_path(self, marker, startnode, endnode, alevel):
		n = startnode
		while n != endnode:
			n.markers[alevel].remove(marker)
			n = n.next[alevel]

	def insert_interval(self, marker, starti, endi):
		if marker in self.lookup:
			raise OverInsertionError()
		startnode = self.insert_node(starti)
		endnode = self.insert_node(endi)
		self.place_marker(marker, startnode, endnode)
		self.lookup[marker] = (starti, endi)

	def place_marker(self, marker, startnode, endnode):
		startnode.starting_markers.add(marker)
		endnode.ending_markers.add(marker)

		starti = startnode.value
		endi = endnode.value

		node = startnode
		alevel = 0

		# mark the non-descending path
		while node.next[alevel].value <= endi:
			# find the level to put the mark on
			while alevel < (node.height - 1) and node.next[alevel+1].value <= endi:
				alevel += 1
			# mark this edge since it's the highest edge that still contains endi
			node.markers[alevel].add(marker)
			node = node.next[alevel]

		# mark non-ascending path
		while node != endnode:
			while alevel > 0 and node.next[alevel].value > endi:
				alevel -= 1
			node.markers[alevel].add(marker)
			node = node.next[alevel]

	def find_containing_node(self, splace):
		markers = set([])
		node = self.head
		for i in range(self.maxheight-1,0,-1):
			# snuggle in
			while node.next[i].value < splace:
				node = node.next[i]
			# drop down
			markers |= node.markers[i]
		# now follow the lowest level 
		while node.next[0].value < splace:
			node = node.next[0]
		markers |= node.markers[0]
		node = node.next[0]
		if node.value == splace:
			markers |= node.starting_markers
		return markers

	def intervals_containing_node(self, splace):
		return [(x, self.lookup[x]) for x in self.find_containing_node(splace)]

	# def find_containing(self, left, right):
	# 	return self.find_containing_node(left) | self.find_containing_node(right)

	def printstuff(self):
		node = self.head
		print 'DUMP'
		while node != self.tail:
			print node.value, '@'+str(node.height), node.starting_markers, node.markers, node.ending_markers
			node = node.next[0]

	def __iter__(self):
		self.current = self.head
		return self

	def next(self):
		while self.current != self.tail and len(self.current.starting_markers)==0:
			self.current = self.current.next[0]
		if self.current == self.tail:
			raise StopIteration
		else:
			me = self.current.starting_markers
			self.current = self.current.next[0]
			return me

	def is_node_contained(self, splace):
		node = self.head
		for i in range(self.maxheight-1,0,-1):
			# snuggle in
			while node.next[i].value < splace:
				node = node.next[i]
			# drop down
			if len(node.markers[i])>0:
				return True
		# now follow the lowest level 
		while node.next[0].value < splace:
			node = node.next[0]
		if len(node.markers[0])>0:
			return True
		node = node.next[0]
		if node.value == splace:
			return len(node.starting_markers)>0
		return False

	# pre: left < right (NOT <=)
	def is_interval_contained(self, left, right):
		node = self.head
		for i in range(self.maxheight-1,0,-1):
			# snuggle in
			while node.next[i].value < left:
				node = node.next[i]
			# drop down
			if len(node.markers[i])>0:
				return True
		# now follow the lowest level 
		while node.next[0].value < left:
			node = node.next[0]
		if len(node.markers[0])>0:
			return True
		node = node.next[0]
		# we are now at the closest point immediately after the left index, which might also be 
		# after the right index as well, so we need to check
		if node.value <= right and len(node.starting_markers)>0:
			return True
		while node.next[0] is not None and node.next[0].value <= right:
			if len(node.starting_markers)>0:
				return True
			node = node.next[0]
		if node.value <= right:
			return len(node.starting_markers)>0
		return False

	def what_contains_interval(self, left, right):
		node = self.head
		markers = set([])
		for i in range(self.maxheight-1,0,-1):
			# snuggle in
			while node.next[i].value < left:
				node = node.next[i]
			# drop down
			markers |= node.markers[i]
		# now follow the lowest level 
		while node.next[0].value < left:
			node = node.next[0]
		markers |= node.markers[i]
		node = node.next[0]
		# we are now at the closest point immediately after or at the left index, which might also be 
		# after the right index as well, so we need to check
		if node.value==left or node.value <= right:
			markers |= node.starting_markers
		while node.next[0] is not None and node.next[0].value <= right:
			markers |= node.starting_markers
			node = node.next[0]
		if node.value <= right:
			markers |= node.starting_markers
		return markers

	def isolated_intervals(self):
		node = self.head
		alive = set([])
		touched = set([])
		isolates = set([])
		while node != self.tail:
			# print node.value, alive, touched, isolates
			if len(node.starting_markers)>0:
				touched |= alive 
			if len(alive)>0:
				touched |= node.starting_markers
			alive |= node.starting_markers
			if len(node.ending_markers)==1 and len(alive)==1 and len(node.ending_markers & touched)==0:
				isolates |= node.ending_markers
			alive -= node.ending_markers
			touched -= node.ending_markers
			node = node.next[0]
		return isolates



"""

execfile('../analysis/skippy.py')
blah = IntervalSkiplist(4)
blah.insert_interval('hi', 0, 5)
blah.insert_interval('hello', 1, 3)
blah.insert_interval('ba', 6, 8)
blah.printstuff()
blah.insert_interval('ab', -3, 0)
blah.isolated_intervals()


execfile('../analysis/skippy.py')
blah = IntervalSkiplist(4)
blah.insert_interval('hi', 0, 5)
blah.insert_interval('hello', 1, 3)
blah.insert_interval('ba', 6, 8)
blah.isolated_intervals()




blah = IntervalSkiplist(4)
blah.insert_interval('hi', 0, 5)
# blah.is_interval_contained(-7,1)
# blah.what_contains_interval(-7,1)
blah.insert_interval('hello', 1, 3)
blah.what_contains_interval(-7,1)
blah.find_containing_node(1)
blah.is_interval_contained(-7,1)


blah = IntervalSkiplist(4)
blah.insert_interval('hi', 0, 5)
print blah.printstuff()
blah.insert_interval('hello', 1, 3)
print blah.printstuff()
print blah.find_containing_node(5)
print blah.find_containing_node(3)
print blah.find_containing_node(-1)
print blah.find_containing(-1,3)

blah.insert_interval('hello', 1, 3)
"""
