function DB = gmm_db_delete_run(db,id_event,id_run)


% Check db for the event
if isfield(db.metadata,'events') && ~intersect(db.metadata.events,id_event)
  warning('Unknown event.')
end