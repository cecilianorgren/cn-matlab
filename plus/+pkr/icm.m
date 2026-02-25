prize_pool = 2.44+1.55+1.03;
prize_pool = 100;
payout = [2.44 1.55 1.03 0 0 0];
payout = payout/sum(payout)*100;
stacks = [4500 3000 1500 725 3775];
nPlayers = numel(stacks);
chance = table(nPlayers,nPlayers);


nIter = 3000;
orig_tickets = [];
order = zeros(nIter,nPlayers);
for kk = 1:nPlayers
    orig_tickets = [orig_tickets; kk*ones(stacks(kk),1)];
end

for pp = 1:nIter
    tickets=orig_tickets;
    for kk=1:nPlayers;
        nTickets = numel(tickets);
        winning_ticket_nr = ceil(nTickets*rand(1,1));
        winner = tickets(winning_ticket_nr);
        order(pp,kk) = winner;
        tickets(tickets==winner) = [];       
    end
end 

probability = zeros(nPlayers,nPlayers);
positions = nPlayers;

for pp = 1:nPlayers;
    [r,c] = find(order==pp);
    for kk = 1:positions
        probability(pp,kk) = numel(find(c==kk));
    end
end
probability=probability/nIter*100;

monetaryValue = zeros(1,nPlayers);
for pp = 1:nPlayers;
    for kk = 1:positions;
        monetaryValue(pp) = monetaryValue(pp)+0.01*probability(pp,kk)*prize_pool*0.01*payout(kk);
    end
end

monetaryValue    