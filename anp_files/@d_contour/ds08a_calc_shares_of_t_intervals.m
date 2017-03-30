function [] = ds08a_calc_shares_of_t_intervals(this)
    
    n_im_poles_distinct = sum([this.im_pz_sorted.type] == 'p');
    n_im_zeros_distinct = sum([this.im_pz_sorted.type] == 'z');
    tf_relative_degree = length(this.poles) - length(this.zeros);
    
    assert(tf_relative_degree >= 0);
    
    weight_crop_inf =       this.weight.of_crop_inf / sqrt(tf_relative_degree+1);
    weight_interpolation =  this.weight.of_axis_parts;
    weight_poles =          this.weight.per_pole * n_im_poles_distinct; % per_pole was 1 before!!
    weight_zeros =          this.weight.per_zero * n_im_zeros_distinct;
    
    weight_sum = weight_crop_inf + weight_interpolation + weight_poles + weight_zeros;
    
    this.shares.crop_inf =  weight_crop_inf / weight_sum;
    this.shares.axis =      weight_interpolation / weight_sum;
    this.shares.poles =  	weight_poles / weight_sum;
    this.shares.zeros =     weight_zeros / weight_sum;    
    
    % now calculate the share between the crop and inf interval
    crop_idx =              strncmp({this.interval_list.type},'crop',4);
    crop_idx =              find(crop_idx == true,1,'first');
    this.arc_lengths.crop = this.interval_list(crop_idx).q_len;
    
    inf_idx =               strncmp({this.interval_list.type},'inf',3);
    inf_idx =               find(inf_idx == true,1,'first');
    this.arc_lengths.inf =  this.interval_list(inf_idx).q_len;
    
    crop_inf_arc_length =   2*this.arc_lengths.crop + this.arc_lengths.inf;
    
    this.shares.crop =      this.arc_lengths.crop/crop_inf_arc_length * this.shares.crop_inf;
    this.shares.inf =       this.arc_lengths.inf /crop_inf_arc_length * this.shares.crop_inf;
end
