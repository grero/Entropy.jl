import Stimulus
import Spiketrains
import JLD

function processTemporalEntropy(sptrains::Dict, bins::Array{Float64,1}, trials::Union(Array{Stimulus.Trial,1},Array{Stimulus.SaccadeTrial,1}),ncombos::Integer, word_size::Integer;skipMissing::Bool=false,alignment::Symbol=:target)
	cells = Spiketrains.sortCells(sptrains)
	tmin = minimum(bins)
	tmax = maximum(bins)
	trial_labels = Stimulus.getTrialLocationLabel(trials,:target)
    D = pmap( _cell -> begin
            fname = string(_cell,"temporalEntropy.jd")
            if isfile(fname)
                println("File already exists. Skipping...")
			else 
                if skipMissing
                    println("Skipping missing data..")
                else
                    isfile(fname) && println("No data found. Computing...")
                    S = sptrains[_cell] 
                    println("Processing cell $_cell")
					spikes = Spiketrains.getTrialRaster(S,trials, alignment,tmin,tmax)
					psth = Spiketrains.getTrialSpikeCount(spikes, bins)
					GTE = Entropy.GroupedTemporalEntropy(int(psth.counts),bins,trial_labels[1:psth.ntrials],word_size)
					JLD.save(fname, {"GTE"=>GTE})
				end
			end
			fname
		end,cells)
	return D
end
