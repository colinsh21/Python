# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:06:30 2015

@author: colinsh
"""

# Imports
import copy
import networkx as nx
import numpy as np
import seaborn; seaborn.set()

class Model(object):
    def __init__(self,graph,initial_infected=1000,reservoir=0,reproduction_number=.65,effective_contacts=15,prob_of_flying_day=.0005,passengers_per_flight=150,screening_effectiveness=0,throttle=0):
        #set pararms
        self.g=graph
        self.initial_infected=initial_infected
        self.reproduction_number=reproduction_number
        self.prob_of_flying_day=prob_of_flying_day
        self.prob_of_flying=1-(1-self.prob_of_flying_day)**30
        self.effective_contacts=effective_contacts
        self.transmission_rate=self.reproduction_number/self.effective_contacts
        self.passengers_per_flight=passengers_per_flight
        self.reservoir=reservoir
        self.screening_effectiveness=screening_effectiveness
        self.throttle=throttle
        
        #set state variables
        self.t=0
        self.history_infected=[]
        self.history_cities_infected=[]
        self.history_infected_interactions=[]
        self.history_space=[]
        self.history_infected_interactions=[]
        self.num_interactions=0
        self.history_num_interactions=[]
        self.time_series_data=[]
        #self call setup
        self.setup_space()
        
    def setup_space(self):
        self.num_infected=0
        self.num_hosts_infected=0
        for node_id in self.g.nodes():
            if node_id == 366: #366 is Jeddah, Saudi Arabia where many MERS infections are
                self.g.node[node_id]['infected']=self.initial_infected
                self.g.node[node_id]["state"] = "I"
                self.num_hosts_infected+=1
                self.num_infected+=self.g.node[node_id]['infected']
                neighbors = self.g.neighbors(node_id)
                for neighbor_id in neighbors: #set throttling for routes
                    print 'before throttle',self.g.edge[node_id][neighbor_id]['weight']
                    self.g.edge[node_id][neighbor_id]['weight']=int(np.ceil(self.g.edge[node_id][neighbor_id]['weight']*(1-self.throttle)))
                    print 'after throttle',self.g.edge[node_id][neighbor_id]['weight']
            else:
                self.g.node[node_id]['infected']=0
                self.g.node[node_id]["state"] = "S"
        
        
        # Track the latest step
        self.history_space.append(copy.deepcopy(self.g))
        self.history_infected.append(self.num_infected)
        self.history_cities_infected.append(self.num_hosts_infected)
        self.num_interactions=0
        self.time_series_data.append((self.num_infected,self.num_hosts_infected,self.num_interactions))
        
        
    def step(self):
        #self.g_step=self.history_space[-1]
        self.g_step=self.g
        self.num_interactions=0
        # Iterate over I and infect any S neighbors
        for self.node_id in self.g.nodes():
            if self.g.node[self.node_id]["state"] == "I":

                # Decide how many infected are traveling
                prob_flying=self.prob_of_flying
                if self.node_id == 366:
                    prob_flying=self.prob_of_flying*(1-self.screening_effectiveness)*(1-self.throttle)
                    
                self.infected_traveling=np.random.binomial(self.g.node[self.node_id]['infected'],prob_flying)

                #check
                if self.infected_traveling>self.g.node[self.node_id]['infected']:
                    self.infected_traveling=self.g.node[self.node_id]['infected']
                    
                # Get travel options
                self.neighbors = self.g.neighbors(self.node_id)
                #print('neighbors',neighbors)
                self.passenger_list=[] #list of passengers on each route
                self.plane_list=[] #list of planes on each route
                for self.neighbor_id in self.neighbors:
                    route_passengers=self.g.edge[self.node_id][self.neighbor_id]['weight']
                    self.passenger_list.append(route_passengers)
                    num_planes=int(np.ceil(float(route_passengers)/self.passengers_per_flight))
                    self.plane_list.append([1.0/float(num_planes)]*num_planes) #probablity of taking a plane
                    #print self.plane_list[-1], 'initialize planes'
                                  

                #print('passenger_list',passenger_list)
                # Create passenger distribution
                self.percent_list = [float(i)/sum(self.passenger_list) for i in self.passenger_list]

                #passenger_dist=np.cumsum(percent_list)
                #print('passenger_dist',percent_list)

                # Create travel distribution
                self.travel_dist=np.random.multinomial(self.infected_traveling, self.percent_list)
                #print('travel_dist',travel_dist)

                # Convert number of travelers to infection
                self.infected_travelers=[]
                for i in xrange(len(self.travel_dist)):
                    self.moving_passengers=self.travel_dist[i]
                    #print self.plane_list[i], 'plane percentage list' 
                    passenger_dist=np.random.multinomial(self.moving_passengers,self.plane_list[i]) #number of infected passengers per plane
                    arriving_infected_total=0
                    #print passenger_dist
                    for j in xrange(len(passenger_dist)):
                         
                        num_infected_plane=passenger_dist[j]
                        p_infect=1-(1-self.transmission_rate)**num_infected_plane #probability a single uninfected passenger gets infected based number of infected on plane
                        arriving_infected_plane=num_infected_plane+np.random.binomial((self.passengers_per_flight-num_infected_plane),
                                                                                      p_infect)
                        #print arriving_infected_plane
                        arriving_infected_total+=arriving_infected_plane
                    
                    self.infected_travelers.append(arriving_infected_total)
                        
                
                    #self.max_passengers=self.g.edge[self.node_id][self.neighbors[i]]['weight']
                    #self.moving_passengers*=self.transmission_rate*self.passengers_per_flight*np.random.normal(1,.15,1)
                    #self.moving_passengers=int(self.moving_passengers)
                    #if self.moving_passengers>self.max_passengers:
                        #self.moving_passengers=self.max_passengers
                        
                    #self.infected_travelers.append(self.moving_passengers)    
                #print self.infected_travelers
                
                # Node update
                # Calculate in host spread
                #print self.node_id,'initial',self.g.node[self.node_id]['infected']
                self.new_infected=int(self.g.node[self.node_id]['infected']*self.reproduction_number)
                #print self.node_id,'spread',self.new_infected
                # Calculate in host death
                self.removed_infected=self.g.node[self.node_id]['infected']
                #print self.node_id,'removed',self.removed_infected
                #print self.node_id,'flying',sum(self.infected_travelers)
                # Spread to neighbors
                for i in xrange(len(self.neighbors)):
                    self.neighbor_id=self.neighbors[i]

                    # Move Infected
                    self.num_infected_at_neighbor=self.g.node[self.neighbor_id]['infected']
                    self.g_step.node[self.neighbor_id]['infected']=self.num_infected_at_neighbor+self.infected_travelers[i]
                    if self.infected_travelers[i]>0:
                        self.history_infected_interactions.append((self.t,self.node_id,self.neighbor_id,self.travel_dist[i],self.infected_travelers[i]))
                        self.num_interactions+=1
                    #self.g_step.node[self.node_id]['infected']=self.g.node[self.node_id]['infected']-self.infected_travelers[i]

                # Sum new infection change
                self.new_infected_total=int(max(self.g.node[self.node_id]['infected']+self.new_infected-self.removed_infected,0))
                
                #maintain reservoir of infection in Jeddah              
                if self.node_id == 366 and self.new_infected_total<self.initial_infected and self.reservoir==1:
                    self.new_infected_total=self.initial_infected
                    #print 'reservoir maintenance'
                    
                self.g_step.node[self.node_id]['infected']=int(self.new_infected_total)
                #print self.node_id,'new infected total',self.new_infected_total
                #print self.node_id,'new infected total graph',self.g_step.node[self.node_id]['infected']
                
        # Update graph
        self.num_infected=0
        self.num_hosts_infected=0
        for self.node_id in self.g_step.nodes():
            if self.g_step.node[self.node_id]['infected']>0:
                #print g.node[self.node_id]['old_label'],'number infected',self.g_step.node[self.node_id]['infected']
                self.g_step.node[self.node_id]["state"] = "I"
                self.num_hosts_infected+=1
                self.num_infected+=int(self.g_step.node[self.node_id]['infected'])
                #print self.node_id,'number infected',self.g_step.node[self.node_id]['infected']
            else:
                self.g_step.node[self.node_id]["state"] = "S"

        # Resolve g_step to g ***Not Resolving***
        self.g=self.g_step

        # Track the latest step
        self.history_space.append(copy.deepcopy(self.g))
        self.t+=1
        self.history_infected.append(self.num_infected)
        self.history_cities_infected.append(self.num_hosts_infected)
        self.history_num_interactions.append(self.num_interactions)
        self.time_series_data.append((self.num_infected,self.num_hosts_infected,self.num_interactions))
        

