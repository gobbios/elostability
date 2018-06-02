# data generation
# need to create sequences that allow defining of rank changes
# ideally:
# - number of rank changes
# - magnitude of rank changes

# start creating a date range on which 'observations were made'
# set inital rank order
# set total observation time per day
# some individuals are more likely to be observed than others
# set steepness of outcome at interaction level (Sanchez-Tojar et al...)
# optional: changes in group composition (maybe leave for later...)



nind <- 7
obsdays <- 15
obsperday <- 4
nchanges <- 2
datagen <- function(nind, obsdays, obsperday, nchanges) {
  # create individual codes
  ids <- sort(apply(expand.grid(letters, letters), 1, paste, collapse = ""))
  ids <- sort(sample(ids, nind))
  # set calender dates
  dates <- seq(from = as.Date("2000-01-01"), by = "day", length.out = obsdays)
  # 'observability' of individuals
  obsb <- runif(n = nind, min = 0.5, max = 1)
  names(obsb) <- ids
  # initial rank orders
  rankmat <- matrix(ncol = nind, nrow = length(dates), 1:nind, byrow = TRUE)
  colnames(rankmat) <- ids
  # select dates for rank changes (and make sure they are neither the first nor the last date)
  changedates <- sort(sample(dates, size = nchanges))
  while(sum(changedates %in% range(dates)) != 0) {
    changedates <- sort(sample(dates, size = nchanges))
  }
  # object to store magnitude of changes
  changemag <- numeric(nchanges)
  # apply rank changes
  if(nchanges >= 1) {
    i = 1
    for(i in 1:nchanges) {
      xdate <- which(dates == changedates[i])
      # who changes
      changers <- sample(ids, 2)
      # ranks before change
      beforechange <- rankmat[xdate - 1, ]
      changerranks <- beforechange[changers]
      # switch
      changedranks <- rev(changerranks)
      # add
      rankmat[xdate, changers] <- changedranks
      # change remainder of rank matrix
      rankmat[xdate:nrow(rankmat), ] <- matrix(ncol = nind, nrow = length(xdate : nrow(rankmat)), rankmat[xdate, ], byrow = TRUE)
      # store magnitude
      changemag[i] <- diff(changerranks)
      # clean up
      rm(xdate, changers, beforechange, changerranks, changedranks)
    }

  }

  # set interactions per day
  iperday <- rpois(n = obsdays, lambda = obsperday)
  mean(iperday)
  if(sum(iperday) < obsdays * obsperday) {
    while(sum(iperday) < obsdays * obsperday) {
      x <- sample(1:length(iperday), 1)
      iperday[x] <- iperday[x] + 1
    }
  }

  if(sum(iperday) > obsdays * obsperday) {
    while(sum(iperday) > obsdays * obsperday) {
      x <- sample(which(iperday > 0), 1)
      iperday[x] <- iperday[x] - 1
    }
  }
  mean(iperday)
  min(iperday)

  # create sequence
  obsdates <- rep(dates, iperday)
  interactants <- matrix(ncol = 2, nrow = length(obsdates), "")
  interactantranks <- matrix(ncol = 2, nrow = length(obsdates), -1)
  for(i in 1:nrow(interactants)) {
    interactants[i, ] <- sample(ids, 2, prob = obsb)
    interactantranks[i, ] <- rankmat[dates == obsdates[i], interactants[i, ]]
  }
  plot(as.numeric(table(interactants)), obsb)


  winnerlosermat <- matrix(ncol = 2, nrow = length(obsdates), "")

  i = 1

  for(i in 1:nrow(winnerlosermat)) {
    maxdiff <- abs(diff(range(rankmat[dates == obsdates[i], ])))
    iids <- interactants[i, ]
    ires <- aniDom:::calculate_winner(rank1 = interactantranks[i, 1], rank2 = interactantranks[i, 2], a = 5, b = 0, max.diff.rank = maxdiff)
    winnerlosermat[i, ] <- interactants[i, ires]
  }

  res <- data.frame(Date = obsdates, winner = winnerlosermat[, 1], loser = winnerlosermat[, 2])

  return(res)
}


datagen(nind = 10, obsdays = 100, obsperday = 4, nchanges = 3)









