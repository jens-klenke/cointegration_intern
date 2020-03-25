#### trendtypes 

# Try Parameter 
trend <- 'none'

if(trend == 'none') {
    ending = -1
    trendype = 1
}else if (trend == 'constant') {
    trendtype = 2
}else if (trend == 'trend') {
    ending = 'FORMULAR' 
    trendtype = 3
}


