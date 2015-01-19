using Winston
import Winston.plot
import Visualizer.visualize

function visualize{T<:String}(::Type{GroupedTemporalEntropy}, fnames::Array{T,1})
	GTE = GroupedTemporalEntropy(fnames[1])
	ta = plot(GTE)

	function func3(_ta,Q,i)
		if i > 0 && i <= length(Q)
			GTE = GroupedTemporalEntropy(Q[i])
			plot(_ta, GTE)
		end
	end
	visualize(fnames,800,600,"Grouped Temporal Entropy", func3, ta)
end

function plot(H::GroupedTemporalEntropy)
	nrows = int(ceil(sqrt(maximum(H.group_labels))))
	ncols = nrows
	ta = Winston.Table(nrows,ncols)
	for i=1:ncols
		for j=1:nrows
		   l = ((i-1)*nrows+j)
		   if (l in H.group_labels) || (i==div(ncols,2)+1) && (j==div(nrows,2)+1)
				ta[j,i] = Winston.FramedPlot()
			end
		end
	end
	ta = plot(ta, H)
	ta
end


function plot(ta::Winston.Table,H::GroupedTemporalEntropy)
	nrows = ta.rows
	ncols = ta.cols
	c = (div(nrows,2)+1 , div(ncols,2)+1)
	ta.attr[:align_interiors] = true
	ta.attr[:cellspacing] = 6.0 #tweaks
	bins = H.bins[1:size(H.H,1)]

	for i=1:ncols
		for j=1:nrows
		   l = ((i-1)*nrows+j)
			if !(l in H.group_labels)
				continue
			end
			idx = find(H.group_labels.==l)
			k = 1
			bidx = find(!isfinite(H.Hc[:,idx]))
			if !isempty(bidx)
				for bi in bidx
					if !(k in bidx)
						Winston.add(ta[j,i],Winston.FillBetween(bins[k:bi-1], H.Hc[:,idx][k:bi-1]-H.dHc[:,idx][k:bi-1],
																bins[k:bi-1], H.Hc[:,idx][k:bi-1] + H.dHc[:,idx][k:bi-1]))
					end
					k = bi + 1

				end
			end
			Winston.add(ta[j,i],Winston.FillBetween(bins[k:end], H.Hc[:,idx][k:end]-H.dHc[:,idx][k:end],
													bins[k:end], H.Hc[:,idx][k:end] + H.dHc[:,idx][k:end]))
														
			Winston.add(ta[j,i],Winston.Curve(bins, H.Hc[:,idx][:]))
		    Winston.setattr(ta[j,i],"title","location $l")
		    Winston.setattr(ta[j,i].x2,"draw_axis",false)
		    Winston.setattr(ta[j,i].y2,"draw_axis",false)
		    Winston.setattr(ta[j,i].x1,"tickdir",1)
		    Winston.setattr(ta[j,i].y1,"tickdir",1)
		end
	end
	#plot the information in the center plot
	p = ta[c...]
	k = 1
	bidx = find(!isfinite(H.I))
	if !isempty(bidx)
		for bi in bidx
			if !(k in bidx)
				Winston.add(p,Winston.FillBetween(bins[k:bi-1], H.I[k:bi-1]-H.dI[k:bi-1],
														bins[k:bi-1], H.I[k:bi-1] + H.dI[k:bi-1]))
			end
			k = bi + 1

		end
	end
	Winston.add(p,Winston.FillBetween(bins[k:end], H.I[k:end]-H.dI[k:end],
											bins[k:end], H.I[k:end] + H.dI[k:end]))
	Winston.add(p, Winston.Curve(bins, H.I))
	Winston.setattr(p.x2,"draw_axis",false)
	Winston.setattr(p.y2,"draw_axis",false)
	Winston.setattr(p.x1,"tickdir",1)
	Winston.setattr(p.y1,"tickdir",1)
	Winston.setattr(p, "xlabel","Time [ms]")
	Winston.setattr(p, "ylabel", "Information [bits]")
	ta
end

function plot(H::TemporalEntropy)
	p = Winston.FramedPlot()
	plot(p, H)
	p
end

function plot(p::Winston.FramedPlot,H::TemporalEntropy)
	Winston.add(p,Winston.Curve(H.bins, H.H))
	Winston.setattr(p.x2,"draw_axis",false)
	Winston.setattr(p.y2,"draw_axis",false)
	Winston.setattr(p.x1,"tickdir",1)
	Winston.setattr(p.y1,"tickdir",1)
	Winston.setattr(p, "xlabel", "Time [ms]")
	Winston.setattr(p, "ylabel", "Entropy [bits]")
	p
end


