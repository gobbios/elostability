# data generation
# need to create sequences that allow defining of rank changes
# ideally:
# - number of rank changes
# - magnitude of rank changes (how exactly? couldn't it be random?)

# start creating a date range on which 'observations were made'
# set inital rank order
# set total observation time per day
# some individuals are more likely to be observed than others
# set steepness of outcome at interaction level (Sanchez-Tojar et al...)
# changes in group composition (migration events)



# nind <- 10
# obsdays <- 15
# obsperday <- 4
# nchanges <- 10
# nmigrations <- 10

datagen <- function(nind, obsdays, obsperday, nchanges, nmigrations = 0) {
  # decide which migrations are imigrations and which are emigrations
  # and a log for migration events
  mig <- NULL
  miglog <- NULL
  if(nmigrations > 0) {
    mig <- sample(c(0, 1), size = nmigrations, replace = TRUE)
    miglog <- data.frame(nb = 1:nmigrations, date = NA, id = NA, mig = c("emig", "imig")[mig+1], imigrank = NA)
  }

  # were there any events?
  events <- nchanges + nmigrations > 0
  if(nchanges + nmigrations > nind) return(NULL)

  # create a log for all events (changes AND migrations)
  if(!events) eventlog <- data.frame(nb = NULL, date = NULL, event = NULL, imigrank = NULL, id1 = NULL, id2 = NULL, magnit = NULL)
  if(events) eventlog <- data.frame(nb = 1:(nchanges + nmigrations), date = NA, event = NA, imigrank = NA, id1 = NA, id2 = NA, magnit = NA)

  # create individual codes (including immigrating individuals)
  ids <- sort(apply(expand.grid(letters, letters), 1, paste, collapse = ""))
  ids <- sort(sample(ids, nind + sum(mig)))
  # set calender dates
  dates <- seq(from = as.Date("2000-01-01"), by = "day", length.out = obsdays)
  # 'observability' of individuals
  obsb <- runif(n = nind + sum(mig), min = 0.33, max = 1)
  obsb <- round(obsb/sum(obsb), 2)
  names(obsb) <- ids
  # initial rank orders (without immigrants)
  rankmat <- matrix(ncol = nind, nrow = length(dates), NA)
  # add imigrants to rank matrix
  if(sum(mig) > 0) {
    im <- which(mig == 1)
    for(i in im) rankmat <- cbind(rankmat, NA)
    rm(im)
  }
  colnames(rankmat) <- ids

  if(events) {
    # select dates for rank changes (and make sure they are neither the first nor the last date)
    if(nchanges > 0) {
      changedates <- sort(sample(dates[2:(length(dates)-1)], size = nchanges))
      eventlog$date[1:length(changedates)] <- sapply(changedates, function(X)which(dates == X))
      eventlog$event[1:length(changedates)] <- "change"
    }

    # select dates for migrations (and make sure they are neither the first nor the last date, nor among the dates of the rank changes)
    if(nmigrations > 0) {
      tempdates <- dates[2:(length(dates) - 1)]
      if(nchanges > 0) tempdates <- tempdates[!tempdates %in% changedates]
      migdates <- sort(sample(tempdates, size = nmigrations))
      miglog$date <- as.character(migdates)
      eventlog$date[(nchanges + 1) : nrow(eventlog)] <- sapply(migdates, function(X)which(dates == X))
      eventlog$event[(nchanges + 1) : nrow(eventlog)] <- c("emig", "imig")[mig + 1]
    }
    # reorder event log
    eventlog <- eventlog[order(eventlog$date), ]
    eventlog$nb <- 1:nrow(eventlog)
    rownames(eventlog) <- NULL
  }

  # presence matrix, which indicates which individuals were present on a given date
  pmat <- matrix(ncol = ncol(rankmat), nrow = length(dates), FALSE)
  colnames(pmat) <- ids
  pmat[, 1:nind] <- TRUE

  # handle migrations in presence matrix
  if(nmigrations >= 1) {
    em <- im <- c()
    # potential emigrators
    if(sum(mig == 0) >= 1) em <- sample(ids[1:nind], size = sum(mig == 0))
    # imigrators
    if(sum(mig) >= 1) im <- ids[(nind +1) : ncol(pmat)]
    # migration dates again
    migdates2 <- migdates
    i = which(eventlog$event != "change")[1]
    for(i in which(eventlog$event != "change")) {
      # temp stuff
      # id <- NULL
      # emigrations
      if(eventlog$event[i] == "emig") {
        pmat[dates >= migdates2[1], em[1]] <- FALSE
        # id <- em[1]
        eventlog$id1[i] <- em[1]
        em <- em[-1]
      }
      # imigrations
      if(eventlog$event[i] == "imig") {
        pmat[dates >= migdates2[1], im[1]] <- TRUE
        # id <- im[1]
        eventlog$id1[i] <- im[1]
        im <- im[-1]
      }
      migdates2 <- migdates2[-1]
      miglog$id[nmigrations - length(migdates2)] <- eventlog$id1[i]
      # rm(id)
    }
    rm(em, im, migdates2)
  }

  # ranks for residents at the beginning
  rankmat[, 1:nind] <- matrix(ncol = nind, nrow = length(dates), 1:nind, byrow = TRUE)
  eventlog

  # now handle events and integrate into rank matrix
  if(events) {
    k = 1
    for(k in 1:nrow(eventlog)) {
      xdate <- eventlog$date[k]
      xrows <- xdate:nrow(rankmat)

      if(eventlog$event[k] == "emig") {
        rankmat[xrows, eventlog$id1[k]] <- NA
        rankmat[xrows, ] <- t(apply(rankmat[xrows, ], 1, rank, na.last = "keep"))
      }
      if(eventlog$event[k] == "imig") {
        xcol <- which(colnames(rankmat) == eventlog$id1[k]) - 1
        # select rank for imigrant:
        imirank <- sample(1:max(rankmat[xdate, ], na.rm = TRUE), 1)
        # miglog$imigrank[i] <- imirank
        eventlog$imigrank[k] <- imirank
        # shift ranks for residents according to where imigrant enters
        oldranks <- na.omit(rankmat[xdate, 1:xcol])
        newranks <- oldranks
        newranks <- newranks[oldranks >= imirank] + 1
        rankmat[xdate, names(newranks)] <- newranks
        rankmat[xdate, xcol + 1] <- imirank
        rankmat[xrows, eventlog$id1[k]] <- imirank
        rankmat[xrows, names(newranks)] <- matrix(ncol = length(newranks), nrow = length(xrows),
                                                  rankmat[xdate, names(newranks)], byrow = TRUE)
        rm(xcol, imirank, oldranks, newranks)
      }
      if(eventlog$event[k] == "change") {
        success <- FALSE
        while(!success) {
          # who changes
          (changers <- sample(ids[pmat[xdate, ]], 2))
          # ranks before change
          beforechange <- rankmat[xdate - 1, ]
          changerranks <- beforechange[changers]
          # check whether changer ranks don't contain NA
          if(!NA %in% changerranks) {
            # switch
            changedranks <- rev(changerranks)
            # add
            rankmat[xdate, changers] <- changedranks
            # change remainder of rank matrix for both individuals
            rankmat[xrows, ] <- matrix(ncol = ncol(rankmat), nrow = length(xrows),
                                              rankmat[xdate, ], byrow = TRUE)
            # fill log
            eventlog$id1[k] <- changers[1]
            eventlog$id2[k] <- changers[2]
            eventlog$magnit[k] <- diff(changerranks)
            # mark as success
            success <- TRUE
          }
        }

        # clean up
        rm(changers, beforechange, changerranks, changedranks, success)

      }

      rm(xdate, xrows)
    }
  }

# interaction generation --------------------------------------------------

  # set interactions per day
  iperday <- rpois(n = obsdays, lambda = obsperday)
  mean(iperday)
  # and make sure the generated data conform to the specified mean
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
    presentids <- pmat[dates == obsdates[i], ]
    presentids <- names(presentids[presentids == TRUE])
    interactants[i, ] <- sample(presentids, 2, prob = obsb[presentids])
    interactantranks[i, ] <- rankmat[dates == obsdates[i], interactants[i, ]]
  }
  # plot(as.numeric(table(interactants)), obsb)


  winnerlosermat <- matrix(ncol = 2, nrow = length(obsdates), "")

  i = 1
  for(i in 1:nrow(winnerlosermat)) {
    maxdiff <- abs(diff(range(rankmat[dates == obsdates[i], ], na.rm = TRUE)))
    iids <- interactants[i, ]
    ires <- aniDom:::calculate_winner(rank1 = interactantranks[i, 1], rank2 = interactantranks[i, 2], a = 0, b = -4, max.diff.rank = maxdiff)
    winnerlosermat[i, ] <- interactants[i, ires]
  }

  res <- data.frame(Date = obsdates, winner = winnerlosermat[, 1], loser = winnerlosermat[, 2])

  # create some more summaries
  xsummary <- data.frame(id = ids, ias = 0, wins = 0, losses = 0, status = "resident", changed = FALSE)
  xsummary$status <- as.character(xsummary$status)
  xsummary$id <- as.character(xsummary$id)
  i = 1
  for(i in 1:nrow(xsummary)) {
    id <- xsummary$id[i]
    xsummary$wins[i] <- sum(winnerlosermat[, 1] == id)
    xsummary$losses[i] <- sum(winnerlosermat[, 2] == id)
    xsummary$ias[i] <- xsummary$losses[i] + xsummary$wins[i]
    if(events) {
      if(id %in% eventlog$id1[eventlog$event == "imig"]) xsummary$status[i] <- "imig"
      if(id %in% eventlog$id1[eventlog$event == "emig"]) xsummary$status[i] <- "emig"
      if(id %in% eventlog$id1[eventlog$event == "change"]) xsummary$changed[i] <- TRUE
      if(id %in% eventlog$id2[eventlog$event == "change"]) xsummary$changed[i] <- TRUE
    }

  }

  reslist <- list(datseq = res, eventlog = eventlog, pmat = pmat, daterange = range(dates), observability = obsb, ids = xsummary)

  if(sum(duplicated(eventlog$date)) > 0) stop("duplicated dates in events...")

  return(reslist)
}


# x <- datagen(nind = 10, obsdays = 20, obsperday = 2, nchanges = 3, nmigrations = 4)
# x$eventlog
# table(x$ids$status, x$ids$changed)
#
# x <- datagen(nind = 10, obsdays = 20, obsperday = 2, nchanges = 0, nmigrations = 4)
# x$eventlog
#
# x <- datagen(nind = 10, obsdays = 20, obsperday = 2, nchanges = 3, nmigrations = 0)
# x$eventlog
#
# x <- datagen(nind = 10, obsdays = 20, obsperday = 2, nchanges = 0, nmigrations = 0)
# x$eventlog
