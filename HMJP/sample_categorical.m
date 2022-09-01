function j = sample_categorical(p)

[p, ind] = sort(p,'descend');

P        = cumsum(p);

j        = ind(find(P(end)*rand<= P,1));
