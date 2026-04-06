add_all_subject_responses_seq = function(des, ys, deads = NULL){
	if (is.null(deads)){
		deads = rep(1, length(ys))
	}
	for (i in seq_along(ys)){
		des$add_one_subject_response(i, ys[i], deads[i])
	}
}
