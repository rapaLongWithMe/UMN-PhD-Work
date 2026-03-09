# How to Clone a Private Repository to Your Local Computer

This guide walks through the process of cloning a private GitHub repository to your local machine. We'll cover both SSH (recommended) and HTTPS methods.

---

## Option 1: SSH (Recommended)

SSH is the preferred method because it doesn't require managing tokens or passwords.

### Step 1: Check for Existing SSH Keys

Open a terminal and check if you already have SSH keys:

```bash
ls -la ~/.ssh
```

You should see files like `id_rsa` and `id_rsa.pub` (or `id_ed25519` and `id_ed25519.pub`).

**If you don't have keys**, generate them:

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

Press Enter for all prompts to use default settings. This creates:
- `~/.ssh/id_ed25519` (private key - keep secret!)
- `~/.ssh/id_ed25519.pub` (public key - safe to share)

### Step 2: Add Your Public Key to GitHub

1. Copy your public key:
   ```bash
   cat ~/.ssh/id_ed25519.pub
   ```

2. Go to GitHub Settings → SSH and GPG keys:
   https://github.com/settings/keys

3. Click "New SSH key"

4. Give it a title (e.g., "MacBook Work Machine")

5. Select key type: **Ed25519** (or RSA if Ed25519 isn't available)

6. Paste the entire public key (starting with `ssh-ed25519` or `ssh-rsa`)

7. Click "Add SSH key"

### Step 3: Test SSH Connection

Verify the connection works:

```bash
ssh -T git@github.com
```

You should see: `Hi [username]! You've successfully authenticated...`

### Step 4: Clone the Repository

Use the SSH URL format (`git@github.com:...`):

```bash
git clone git@github.com:username/repo-name.git ~/path/to/destination
```

**Example:**
```bash
git clone git@github.com:rapaLongWithMe/UMN-PhD-Work.git ~/Desktop/Git/UMN-PhD-Work
```

---

## Option 2: HTTPS with Personal Access Token

Use this method if you prefer HTTPS over SSH.

### Step 1: Generate a Personal Access Token on GitHub

1. Go to GitHub Settings → Developer settings → Personal access tokens:
   https://github.com/settings/tokens/new

2. Give it a name (e.g., "Local Machine Access")

3. Set expiration (recommend 90 days or custom)

4. Select scopes: check the `repo` box (full control of private repositories)

5. Click "Generate token"

6. **Copy the token immediately** - you won't see it again!

### Step 2: Set Up Git Credential Helper

On macOS, use the keychain to store credentials:

```bash
git config --global credential.helper osxkeychain
```

### Step 3: Clone the Repository

Use the HTTPS URL format:

```bash
git clone https://github.com/username/repo-name.git ~/path/to/destination
```

**Example:**
```bash
git clone https://github.com/rapaLongWithMe/UMN-PhD-Work.git ~/Desktop/Git/UMN-PhD-Work
```

When prompted:
- Username: your GitHub username
- Password: paste your personal access token (not your GitHub password!)

The credentials will be saved to your keychain for future use.

---

## Common Issues and Solutions

### Issue: "Repository not found"

**Cause:** Wrong URL, repository doesn't exist, or no access

**Solutions:**
- Double-check the repository URL
- Verify the repository is actually private and you have access
- For private repos: confirm you're using the correct authentication method

### Issue: "Permission denied (publickey)" with SSH

**Cause:** SSH key not authorized on GitHub

**Solutions:**
1. Verify your public key is added to GitHub:
   ```bash
   ssh -T git@github.com
   ```

2. Check that SSH is using the correct key:
   ```bash
   ssh-add -l
   ```

3. If your key isn't listed, add it:
   ```bash
   ssh-add ~/.ssh/id_ed25519
   ```

4. Try the connection again

### Issue: HTTPS authentication fails repeatedly

**Cause:** Token expired, wrong credentials, or keychain issues

**Solutions:**
1. Generate a new personal access token
2. Clear old credentials from keychain:
   ```bash
   git credential-osxkeychain erase
   # Type: host=github.com
   # Press Ctrl+D twice
   ```
3. Try cloning again to re-enter credentials

### Issue: "Key is invalid. You must supply a key in OpenSSH public key format"

**Cause:** Wrong key format pasted into GitHub

**Solutions:**
- Make sure you're pasting the **public key** (from `.pub` file), not the private key
- Ensure the entire key is copied with no extra spaces or line breaks
- Key should start with `ssh-ed25519` or `ssh-rsa`

---

## After Cloning: Verifying Your Setup

Once cloned, verify everything is working:

```bash
cd ~/path/to/destination
git status                    # Check branch and status
git remote -v                 # Verify remote URL
git log --oneline -5          # View recent commits
```

---

## Useful Git Commands

```bash
# Check current remote
git remote -v

# Change remote URL (if needed)
git remote set-url origin git@github.com:username/repo.git

# Verify authentication
ssh -T git@github.com

# Configure git user (if needed)
git config user.name "Your Name"
git config user.email "your@email.com"

# Make these settings global (apply to all repos)
git config --global user.name "Your Name"
git config --global user.email "your@email.com"
```

---

## Summary: Quick Checklist

- [ ] SSH keys exist or generated (`~/.ssh/id_ed25519.pub`)
- [ ] Public key added to GitHub (Settings → SSH and GPG keys)
- [ ] SSH connection tested (`ssh -T git@github.com`)
- [ ] Repository URL copied (use SSH format: `git@github.com:...`)
- [ ] Clone command executed
- [ ] Repository cloned successfully to desired location
- [ ] Verified with `git status` and `git remote -v`

