function combineTwoChannelStacks_tile(path1, path2, outPath)
    m1 = Movie.load(path1);
    m2 = Movie.load(path2);
    assert(m1.nX==m2.nX && m1.nY==m2.nY && m1.nFrames==m2.nFrames);
    n = m1.nFrames;
    for k = 1:n
        f1 = m1.readFrames(k);
        f2 = m2.readFrames(k);
        tiled = cat(2, f1, f2);
        if k == 1
            imwrite(tiled, outPath, 'Compression', 'none');
        else
            imwrite(tiled, outPath, 'WriteMode', 'append', 'Compression', 'none');
        end
        if mod(k, 500) == 0, fprintf('frame %d / %d\n', k, n); end
    end
    fprintf('Wrote %s (%d frames, %d x %d)\n', outPath, n, size(tiled,1), size(tiled,2));
end