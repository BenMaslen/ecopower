utils::globalVariables(c(
  "object",
  "N",
  "coeffs",
  "newdata",
  "n_replicate",
  "do.fit",
  "test",
  "term",
  "object0"
))

MVApowerstat = function(stats, coeffs) {
  anova(
    extend(
      object=object,
      N=N,
      coeffs=coeffs,
      newdata=newdata,
      n_replicate=n_replicate,
      do.fit=do.fit
    ),
    nBoot=1,
    test=test,
    show.time = "none"
  )$table[term,3]
}

MVApowerstatObj = function(stats) {
  object_sim = extend(
    object=object,
    N=N,
    coeffs=coeffs,
    newdata=newdata,
    n_replicate=n_replicate,
    do.fit=do.fit
  )

  object0_sim = update(
    object_sim,
    formula=object0$formula[-2],
    data=object_sim$data
  )

  anova(
    object0_sim, object_sim, nBoot=1, test=test, show.time="none"
  )$table[2,3]
}
