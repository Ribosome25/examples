from mpi4py import MPI
import numpy as np
import pandas as pd
import pickle

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
n_processors = comm.Get_size()
if my_rank == 0:
    print("Processors found: ", n_processors)

def scatter_list_to_processors(comm, data_list, n_processors):
    import math
    data_amount = len(data_list)
    heap_size = math.ceil(data_amount/(n_processors-1))

    for pidx in range(1, n_processors):
        try:
            heap = data_list[heap_size*(pidx-1):heap_size*pidx]
        except IndexError:
            heap = data_list[heap_size*(pidx-1):]
        comm.send(heap, dest=pidx)
    return True

def receive_from_processors_to_list(comm, n_processors):
    # receives list, combine them and return
    feedback = []
    for pidx in range(1,n_processors):
        received = comm.recv(source=pidx)
        feedback.append(received)
    return feedback

  def receive_from_processors_to_dict(comm, n_processors):
    # receives dicts, combine them and return
    feedback = dict()
    for pidx in range(1,n_processors):
        received = comm.recv(source=pidx)
        feedback.update(received)
    return feedback
  
  if my_rank == 0:
    for i_iter in iters:
        # Update something to all processors. 
        broadcast_msg = update 
        comm.bcast(broadcast_msg, root=0)
        # send task list to all processors. 
        scatter_list_to_processors(comm, task_list, n_processors)
        result_dict = receive_from_processors_to_dict(comm, n_processors)
        do some thing... 
  else:
    for i_iter in iters:
        broadcast_msg = update    # claim a size

        received_update = comm.bcast(broadcast_msg, root=0)
        task_list = comm.recv(source=0)
        do something...
        comm.send(results, dest = 0)
        
