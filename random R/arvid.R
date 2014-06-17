getWindowValues <- function (data, window = 50000, breakSize = 1000){
	start = data$Crub_start[1]
	pCend = data$Crub_end[1]
	end = start + window
	
	myX = c()
	rY = c()
	pcEnd = data$Crub_end[1]
	#Rlength = pCend-start
	
	tY = c()
	pTend = data$Ath_end[1]
	aTlength = data$Ath_end[1]-data$Ath_start[1]
	
	
	lY = c()
	pLend = data$Alyr_end[1]
	aLlength = data$Alyr_end[1]-data$Alyr_start[1]
	
	while (i < length(data$Crub_start)){
		if (data$Crub_end[i] > end){
			lY = append(lY, aLlength)
			tY = append(tY, aTlength)
			myX = append(myX, (start+pCend)/2)
			rY = append(rY, pCend-start)
			aLlength = 0
			aTlength = 0
			Rlength = 0
			start = data$Crub_start[i]
			end = start + window
			pCend = data$Crub_end[i]
		} else{
			aLlength = aLlength + (data$Alyr_end[i]-data$Alyr_start[i])
			aTlength = aTlength + (data$Ath_end[i]-data$Ath_start[i])
			Rlength = Rlength + (data$Crub_end[i]-data$Crub_start[i])
			
			tmpL = data$Alyr_start[i]-pLend
			tmpT = data$Ath_start[i]-pTend
			
			if (!(tmpL < 0 | tmpL > breakSize)){
				aLlength = aLlength + tmpL
				#Rlength = Rlength + (pcEnd-data$Crub_start[i])
			} 
			if (!(tmpT < 0 | tmpT > breakSize)){
				aLlength = aLlength + tmpT
				#Rlength = Rlength + (pcEnd-data$Crub_start[i])
			}
			
			pcEnd = data$Crub_end[i]
			pLend = data$Alyr_end[i]
			pTend = data$Ath_end[i]
		}
		
		i = i + 1
	}
	
	list(myX, tY, lY, rY)
}