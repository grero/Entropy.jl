import NSBEntropy
import Spiketrains
import Information
import StatsBase

function example1()
    nreps = 1000
    #20 ms window
    t = linspace(0,.02,1000)
    dt = mean(diff(t))
    #constant firing rate of 260 Hz
    fr = 260.*ones(length(t))
    T,repidx = Spiketrains.generateSpikes(fr,dt,nreps)
    #create an aligned spiketrain
    at1 = Information.AlignedSpiketrain(T,repidx,nreps,length(T))
    #create a histogram using 0.5 ms bins
    countsT = Information.getTrialSpikeCount(at1,[0:0.0005:0.015])
    #ensure that we have a maxiumum of 1 spike per bin
    countsT.counts[countsT.counts.>1] = 1
    #hask the words
    Q = countsT.counts'*(2.^[0:29])
    #get the counts of each unique value
    q = StatsBase.countmap(Q)
    #compute entropy, making use of the fact that the refractory period ensures L <= 10^16
    Entropy.NSBEntropy.saddlepoint(collect(values(q)),10^16)
end
