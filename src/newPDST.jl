durations = [1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 7.5, 9.0, 10.5, 12.0, 12.5, 15.0]
min_periods = [0.5, 0.52, 0.6517, 0.7824, 0.912, 1.178, 1.3056, 1.567, 1.952, 2.343, 2.75, 3.14, 3.257, 3.91]
max_periods = [50.045, 118.626, 231.69, 400.359, 635.76, 725, 725, 725, 725, 725, 725, 725, 725, 725]
num_dur = 14

function get_legal_durations(period::Float64,duration::Float64)
  min_duration = 0.0
  max_duration = 0.0
  i = 1
  #determine what maximum and minimum durations were searched for this period
  while min_duration == 0.0 || max_duration == 0.0
    if i > 14
      println("No durations match this period")
      return(0.0,0.0)
    end    
    if period <= max_periods[i] && min_duration == 0.0
      min_duration = durations[i]
    end 
    if period >= min_periods[num_dur+1-i] && max_duration == 0.0
      max_duration = durations[num_dur+1-i]
    end
    i+=1
  end
  if duration<=max_duration/24 && duration>=min_duration/24
    return duration
  elseif duration>=max_duration/24
    return max_duration/24
  elseif duration <=min_duration/24
    return min_duration/24
  end
end
