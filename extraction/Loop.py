import itertools
from extraction.Experiment import PhotonType, Photon
from time import sleep

import random, math
import perceval as pcvl
import perceval.components as symb
import perceval.utils as pcvlutils
import multiprocessing as mp
import re, signal

class Loop(object):

    def __init__(self, width, depth, phs_in, phs_out, loss_Source=0, loss_loop_inner=0, loss_loop_outer=0, in_state=[],
                 indistinguishability=1, reflectivity_pattern=1):
        self.loss_Source=loss_Source
        self.loss_loop_inner=loss_loop_inner
        self.loss_loop_outer=loss_loop_outer

        self.phs_in = phs_in
        self.phs_out = phs_out
        self.width=width
        self.depth = depth

        self.reflectivity_pattern = reflectivity_pattern

        self.m_in = width
        self.m_out = width
        self.in_state = None
        self.out_states = {}

        self.experiment_circuit = None

        self.beamsplitter = []
        self.photons = {}
        Photon.newid = itertools.count()

        self.m_loss = self.interferometer(width, depth)

        for i in range(phs_in):
            self.photons[i] = Photon(type=PhotonType.READOUT)
        for i in range(self.m_loss):
            self.photons[phs_in + i] = Photon(type=PhotonType.LOSS)

        if in_state != []:
            self.calc_in_and_out_states_with_losses(in_state=in_state,indistinguishability=indistinguishability)

    def calculate_input_state_with_losses(self, in_state, indistinguishability=1):
        if indistinguishability == 0:
            in_state = in_state + [0] * self.m_loss
            self.in_state = '|' + ",".join([f'{in_state[i]}{f"{{_:{i}}}" if in_state[i] > 0 else ""}' for i in range(len(in_state))]) + '>'
        elif indistinguishability == 1:
            self.in_state = "|" + ",".join([f'{s}' for s in in_state + [0] * self.m_loss]) + ">"
        return self.in_state

    def calc_out_states(self, fore_calc=False):
        assert self.in_state != None, "In State is required"

        if not fore_calc and self.out_states != {}:
            return self.out_states

        postprocess = lambda o: (sum(o[i] > 0 for i in range(self.m_out)) == self.phs_out
                                 and sum(o[i] == 1 for i in range(self.m_out, self.m_out+self.m_loss)) <= self.phs_in - self.phs_out
                                 and sum(o) == self.phs_in
                                 and (self.m_loss == 0 or self.phs_in - self.phs_out == 0 or max(o[i] for i in range(self.m_out,self.m_out+ self.m_loss)) <= 1)
                                 )

        out_states = {}
        for st in itertools.combinations(range(self.m_out), self.phs_out):
            out = [0] * self.m_out + [0] * self.m_loss
            for i in st:
                out[i] = 1
            #state = pcvl.BasicState(out + [0] * self.m_loss)
            st = "|" + ",".join([f'{s}' for s in out]) + ">"
            if postprocess(out):
                out_states[tuple(out)] = st
            if self.phs_in - self.phs_out > 0:
                losses = partitions(self.phs_in - self.phs_out)
                for l in losses:
                    for st2 in itertools.combinations(range(self.m_out + self.m_loss), len(l)):
                        out2 = out.copy()
                        j = 0
                        for i2 in st2:
                            out2[i2]+=l[j]
                            j+=1
                        if postprocess(out2):
                            out_states[tuple(out2)] = st
        self.out_states = out_states
        return self.out_states


    def calc_in_and_out_states_with_losses(self,in_state=[], indistinguishability=0):
        if in_state == [] or in_state == None:
            in_state = [x for x in pcvl.BasicState(self.in_state)][:self.width]
        self.calculate_input_state_with_losses(in_state, indistinguishability)
        self.calc_out_states()

    def calculate_no_wires(self, width, depth):
        m = width
        if self.loss_Source != 0 and self.loss_loop_outer != 0:
            if self.loss_loop_inner != 0:
                m += (width -1) * depth #outer loop + source
                m += width * depth # inner loops
            else:
                m += width * depth
        elif self.loss_Source != 0:
            if self.loss_loop_inner != 0:
                m += width - 1 #source
                m += width * depth # inner loops
            else:
                m += width #source
        elif self.loss_loop_outer != 0:
            if self.loss_loop_inner != 0:
                m += (depth-1) * (width - 1) # outer loops
                m += width * depth # inner loops
            else:
                m += (depth-1) * (width)
        elif self.loss_loop_inner != 0:
            m += width * depth  # inner loops
        return m

    def reflectivity(self,m):
        if self.reflectivity_pattern==2 and self.m_in == 10: #special lorenzo pattern
            return list([.5,.6,.7,.8,.8,.8,.7,.5,.4])[m]
        else:
            return (m + 2 - 1) / (m + 2)

    def interferometer(self,width, depth):
        m = self.calculate_no_wires(width, depth)

        circ = symb.Circuit(m=m, name="Interferometer")
        k = 0
        for i in range(depth):
            for j in range(width):
                if i == 0 and j == 0:  # first photon in first loop
                    k += addLoss(circ, j, width + k, self.loss_Source + self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_Source)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                elif i == 0 and j > 0 and j < width - 2:  # second+ photon in first loop
                    k += addLoss(circ, j, width + k, self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_Source)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                elif i == 0 and j > 0 and j == width - 2:  # last photon in first loop
                    k += addLoss(circ, j, width + k, self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_Source)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                    k += addLoss(circ, j + 1, width + k, self.loss_loop_inner)
                elif i > 0 and j == 0: # first photons in second+ loop
                    k += addLoss(circ, j, width + k, self.loss_loop_outer + self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_loop_outer)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                elif i > 0 and j < width - 2:  # last photon in second+ loop
                    k += addLoss(circ, j, width + k, self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_loop_outer)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                elif i > 0 and j == width - 2:  # all photons in second+ loop
                    k += addLoss(circ, j, width + k, self.loss_loop_inner)
                    k += addLoss(circ, j + 1, width + k, self.loss_loop_outer)
                    circ.add(j, symb.BS(symb.BS.r_to_theta(self.reflectivity(j))))
                    k += addLoss(circ, j + 1, width + k, self.loss_loop_inner)
            circ.add(0, symb.PERM(list(range(width))))
        self.experiment_circuit = circ
        return k

    def rowNormEstimator(self, experimental_output, in_state = None):

        assert self.in_state != None or in_state != None, "In State is required"
        assert experimental_output != {}, "Experimental output should be a dictionary"

        if in_state != None:
            self.calculate_input_state_with_losses(in_state=in_state, indistinguishability=1)

        in_state = pcvl.BasicState(self.in_state)

        u = self.experiment_circuit.compute_unitary()
        shuffled_out = []
        for k,v in experimental_output.items():
            for _ in range(v):
                shuffled_out.append(k)
        random.shuffle(shuffled_out)

        p = lambda ins, outs: math.prod([sum([abs(u[l][k])**2 for k in range(len(outs)) if outs[k] != 0]) for l in range(len(ins)) if ins[l] != 0])

        out = [0]
        for k in shuffled_out:
            p_k = p(in_state,k)
            out.append(out[len(out)-1]+1 if p_k >= (self.phs_in/(self.m_out+self.m_loss))**self.phs_in else out[len(out)-1]-1)

        return out



    def run(self,in_state = [], processes = 1, chunks = 1, backend="Naive", indistinguishability=1, t_sleep=10):

        assert self.in_state != None or self.in_state != [] or in_state != None or in_state != [], "In State is required"

        self.calc_in_and_out_states_with_losses(in_state, indistinguishability)
        dist = {i:0 for i in self.out_states.values()}
        if processes > 1:
            with mp.Manager() as manager:
                distribution = manager.dict({i:0 for i in set(self.out_states.values())})
                counter = manager.Value('i', 0)
                shutdown = manager.Event()
                signal.signal(signal.SIGTERM, lambda signum, frame: shutdown.set())

                lock = manager.Lock()
                object_args = [self.width, self.depth, self.phs_in, self.phs_out, self.loss_Source, self.loss_loop_inner, self.loss_loop_outer]
                args = [(self.in_state, self.out_states, indistinguishability, distribution, counter, backend, object_args, chunks, lock, i, shutdown) for i in range(chunks)]
                with manager.Pool(processes=processes) as pool:
                    pool.map_async(Loop.sim, args)
                    while True:
                        c = 0
                        with lock:
                            c = counter.value
                        print("calculated {} of {}".format(c, len(self.out_states)), flush=True)
                        if c == len(self.out_states):
                            break
                        sleep(t_sleep)
                dist=dict(distribution)
        else:
            backend = pcvl.BackendFactory.get_backend(backend)
            sim = backend(self.experiment_circuit)
            i=0
            for k, v in self.out_states.items():
                #print("calculate {} of {}".format(i,len(self.out_states)))
                dist[v] += sim.prob(input_state=pcvl.BasicState(self.in_state), output_state=pcvl.BasicState(k))
                i+=1
        return dist

    @staticmethod
    def sim(args) -> None:
        (in_state, out_states, indistinguishability, dist, counter, backend, loopArgs, workers, lock, id, shutdown) = args
        loop = Loop(*loopArgs)
        loop.in_state = in_state
        loop.out_states = out_states

        backend = pcvl.BackendFactory.get_backend(backend)
        sim = backend(loop.experiment_circuit)
        (first, last) = (int(len(loop.out_states)/workers*id),int(len(loop.out_states)/workers*(id+1)))
        out = dict(itertools.islice(loop.out_states.items(),first,last))
        for k,v in out.items():
            p = sim.prob(input_state=pcvl.BasicState(loop.in_state), output_state=pcvl.BasicState(k))
            #print("Worker {} update for {}".format(id, pcvl.BasicState(k)), flush=True)
            with lock:
                dist[v] += p
                counter.value += 1
                if shutdown.is_set():
                    #print("Have to shutdown", flush=True)
                    return

def addLoss(circ,m1,m2,r,merge=False):
    if r>0:
        # circ.add_modes(m2-1,1)
        circ.add(m1,loss(m1,m2,r=r),merge=merge)
        return 1
    return 0

def loss(m1,m2,r):
    m = abs(m2-m1)
    loss = symb.Circuit(m+1, "Loss m1/m2={}/{} r={}".format(m1,m2,1-r))
    if m > 1:
        a1,a2,*b,c = list(range(0,m+1))
        perm=symb.PERM([a1,c,*b,a2])
        loss.add(0,perm)
    loss.add(0, symb.BS(symb.BS.r_to_theta(1-r)))
    if m > 1:
        loss.add(0,perm)
    return loss

def lossyBS(r,l_in_1=0,l_out_1=0,l_in_2=0,l_out_2=0):
    BS = pcvl.Processor("Naive",pcvl.Circuit(2,name="Loop BS"))
    if l_in_1 != 0: BS.add(0,pcvl.LC(l_in_1))
    if l_in_2 != 0: BS.add(1,pcvl.LC(l_in_2))
    BS.add(0,pcvl.BS(pcvl.BS.r_to_theta(r)))
    if l_out_1 != 0: BS.add(0,pcvl.LC(l_out_1))
    if l_out_2 != 0: BS.add(1,pcvl.LC(l_out_2))
    return BS

def calc_out_states(in_phs, out_phs, m_out, m_lcs):
    out_states = {}
    for st in itertools.combinations(range(m_out), in_phs):
        l = [0] * m_out + [0] * m_lcs
        l2 = [0] * m_out + [0] * m_lcs
        for i in st:
            l[i] = 1
            l2[i] = '?'
        state = pcvl.BasicState(l)
        out_states[state] = '|' + ",".join([f'{ph}' for ph in l2]) + '>'
    return out_states

def partitions(n, I=1):
    yield (n,)
    for i in range(I, n//2 + 1):
        for p in partitions(n-i, i):
            yield (i,) + p

